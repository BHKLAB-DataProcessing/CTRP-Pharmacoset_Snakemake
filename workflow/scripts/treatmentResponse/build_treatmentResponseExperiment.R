## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

    save.image(paste0(snakemake@rule, ".RData"))
}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(CoreGx))

####################################################################################
# 0.2 read in metadata
# --------------------

zipDir <- tools::file_path_sans_ext(INPUT$tr)
dir.create(zipDir, recursive = TRUE, showWarnings = FALSE)
print(paste0("Unzipping ", INPUT$tr, " into ", zipDir))
unzip(INPUT$tr, exdir = zipDir)

files <- paste0(zipDir, "/", list.files(zipDir))


# Get the treatment metadata
# --------------------------
treatmentMetadata <- 
    files[grepl(pattern = "*per_compound*", x = files)] |>
    fread() 
treatmentMetadata[, treatmentid := cpd_name]

# Get the sample metadata
# -----------------------
sampleMetadata <- 
    files[grepl(pattern = "*per_cell_line*", x = files)] |>
    fread()

sampleMetadata <- sampleMetadata[, .(
        sampleid = ccl_name, 
        tissueid = ccle_primary_site, 
        master_ccl_id = as.character(master_ccl_id), 
        sample_availability = ccl_availability, 
        ccle_primary_hist = ccle_primary_hist)]


# Read in sensitivity data
# ------------------------
sensitivityRaw <- 
    files[grepl(pattern = "*per_cpd_post_qc*", x = files)] |>
    fread()


# code taken from https://github.com/BHKLAB-DataProcessing/PSet_CTRPv2-snakemake/blob/master/Oldscripts/downloadSensData.R
ctrp.sensitivityInfo <-  read.delim(files[grepl(pattern = "*meta.per_experiment*", x = files)])

repExps <- unique(ctrp.sensitivityInfo$experiment_id[duplicated(ctrp.sensitivityInfo$experiment_id)])

## Just collapsing the dates, everything else is identical.
for (exp in repExps) {
  myx <- ctrp.sensitivityInfo$experiment_id == exp
  duplicates <- duplicated(ctrp.sensitivityInfo$experiment_id) & myx
  first <- myx & !duplicates
  # print(ctrp.sensitivityInfo[myx,])
  # browser()
  ctrp.sensitivityInfo[first, ] <- apply(ctrp.sensitivityInfo[myx, ], 2, function(x) paste(unique(x), collapse = "//"))

  ctrp.sensitivityInfo <- ctrp.sensitivityInfo[!duplicates, ]
}


sensitivityInfo <- data.table::as.data.table(ctrp.sensitivityInfo)


sensitivityInfo <- merge(
    sensitivityInfo,
    sampleMetadata,
    by = "master_ccl_id",
    all.x = TRUE
)


sensitivityRaw <- merge(
    sensitivityRaw,
    sensitivityInfo[, .(sampleid, experiment_id = as.integer(experiment_id), culture_media)],
    by = "experiment_id"
)

sensitivityRaw <- merge(
    sensitivityRaw,
    treatmentMetadata[, .(treatmentid, master_cpd_id)],
    by = "master_cpd_id"
)

rawdt <- sensitivityRaw[, .(
    sampleid, culture_media, 
    treatmentid, experiment_id, 
    dose = cpd_conc_umol, viability =  100*cpd_avg_pv
)]

# Add technical replicate column 
# rawdt[, .N, by = .(sampleid, treatmentid, culture_media, dose)][order(N)]

rawdt[,
    tech_rep := seq_len(.N),
    by=.(sampleid, treatmentid, culture_media, dose)
]



# -- Create a TREDataMapper
print("Creating a TREDataMapper")
treDataMapper <- CoreGx::TREDataMapper(rawdata = rawdt)

# groups <- list(
#     rowDataMap = c("treatmentid", "dose"),
#     rowDataMapTry2 = c("treatmentid", "dose", "tech_rep", "experiment_id"),
#     colDataMap = c("sampleid", "culture_media"),
#     assayMap1 = c("treatmentid", "dose", "sampleid", "culture_media"),
#     assayMap = c("treatmentid", "dose", "tech_rep", "experiment_id", "sampleid", "culture_media")
# )


# subsets <- list(FALSE, TRUE, FALSE, FALSE, TRUE)
# guess <- guessMapping(treDataMapper, groups=groups, subset=FALSE)
# guess
# # guess didnt work so well... 

CoreGx::rowDataMap(treDataMapper) <- list(
    id_columns = c("treatmentid", "dose", "tech_rep"),
    mapped_columns = c()
)

CoreGx::colDataMap(treDataMapper) <- list(
    id_columns = c("sampleid", "culture_media"),
    mapped_columns = c()
)

CoreGx::assayMap(treDataMapper) <- list(
    sensitivity = list(
        id_columns = c("treatmentid", "dose", "tech_rep", "sampleid", "culture_media"),
        mapped_columns = c("viability")
    )
)

print("Running CoreGx::metaConstruct")
(tre <- metaConstruct(treDataMapper))

print("Endoaggregating to create mono_viability Assay")
tre_qc <- tre |>
    endoaggregate(
        assay="sensitivity",
        target="sensitivity",  # create a new assay named mono_viability
        mean_viability=pmin(100, mean(viability)), # pmin takes the minimum of two vectors element-wise
        by=c("treatmentid", "dose", "sampleid"),
        nthread=THREADS
)

print("Endoaggregating to create profiles_recomputed Assay")
tre_fit <- tre_qc |> CoreGx::endoaggregate(
    {  # the entire code block is evaluated for each group in our group by
        # 1. fit a log logistic curve over the dose range
        fit <- PharmacoGx::logLogisticRegression(dose, mean_viability,
            viability_as_pct=TRUE)
        # 2. compute curve summary metrics
        ic50 <- PharmacoGx::computeIC50(dose, Hill_fit=fit)
        aac <- PharmacoGx::computeAUC(dose, Hill_fit=fit)
        # 3. assemble the results into a list, each item will become a
        #   column in the target assay.
        list(
            HS=fit[["HS"]],
            E_inf = fit[["E_inf"]],
            EC50 = fit[["EC50"]],
            Rsq=as.numeric(unlist(attributes(fit))),
            aac_recomputed=aac,
            ic50_recomputed=ic50
        )
    },
    assay="sensitivity",
    target="profiles",
    enlist=FALSE,  # this option enables the use of a code block for aggregation
    by=c("treatmentid", "sampleid"),
    nthread=THREADS  # parallelize over multiple cores to speed up the computation
)

print(paste("Saving Output Files to", OUTPUT$tre))

saveRDS(tre_fit, file = OUTPUT$tre)
saveRDS(sensitivityInfo, file = OUTPUT$sensInfo)