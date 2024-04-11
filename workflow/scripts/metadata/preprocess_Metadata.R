## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads

    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(
            file = snakemake@log[[1]], 
            append = FALSE, 
            type = c("output", "message"), 
            split = TRUE
    )

    file.path("resources", paste0(snakemake@rule, ".RData")) |> 
        save.image()
}else{
    file.path("resources", "preprocess_Metadata.RData") |>
        load()
}




#############################################################################
# Load INPUT
#############################################################################

zipDir <- file.path(snakemake@resources$tmpdir, snakemake@rule, tools::file_path_sans_ext(INPUT$tr))
dir.create(zipDir, recursive = TRUE, showWarnings = FALSE)
print(paste0("Unzipping ", INPUT$tr, " into ", zipDir))
unzip(INPUT$tr, exdir = zipDir)
files <- paste0(zipDir, "/", list.files(zipDir))


#############################################################################
# Main Script
# Simply load in the metadata files and then save them as tsv files

message("Reading in the metadata files")

# Get the treatment metadata
# --------------------------
message("Reading in treatment metadata")
treatmentMetadata <- 
    files[grepl(pattern = "*per_compound*", x = files)] |>
    data.table::fread() 
treatmentMetadata[, treatmentid := cpd_name]
str(treatmentMetadata)


# Get the sample metadata
# _______________________
# v20.meta.per_cell_line
message("Reading in sample metadata")
sampleMetadata <- 
    files[grepl(pattern = "*per_cell_line*", x = files)] |>
    data.table::fread()

sampleMetadata <- sampleMetadata[, .(
        sampleid = ccl_name, 
        tissueid = ccle_primary_site, 
        master_ccl_id = as.character(master_ccl_id), 
        sample_availability = ccl_availability, 
        ccle_primary_hist = ccle_primary_hist)]

str(sampleMetadata)

## --------------------- Save OUTPUT ------------------- ##

message("Saving the metadata files")
message("Saving treatment metadata to ", OUTPUT$rawTreatmentMetadata)

dir.create(dirname(OUTPUT$rawTreatmentMetadata), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(
    treatmentMetadata,
    file.path(OUTPUT$rawTreatmentMetadata),
    sep = "\t",
    quote = FALSE
)

message("Saving sample metadata to ", OUTPUT$rawSampleMetadata)
dir.create(dirname(OUTPUT$rawSampleMetadata), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(
    sampleMetadata,
    file.path(OUTPUT$rawSampleMetadata),
    sep = "\t",
    quote = FALSE
)
