## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    RESOURCES <- snakemake@resources    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

    save.image(paste0(snakemake@rule, ".RData"))
}


zipDir <- file.path(RESOURCES$tmpdir, snakemake@rule, tools::file_path_sans_ext(INPUT$tr))
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

ann_cols <- treatmentMetadata[, .(cpd_name, broad_cpd_id)]


(compound_nameToCIDS <- AnnotationGx::getPubchemCompound(
    ann_cols[, cpd_name],
    from='name',
    to='cids',
))

properties=c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES')
message(
    "Getting the following properties from PubChem: ", 
    paste(properties, collapse= " "), " for ", nrow(compound_nameToCIDS), " compounds")

pubchemProperties <- 
    compound_nameToCIDS[
        1:nrow(compound_nameToCIDS), 
        AnnotationGx::getPubchemCompound(ids = cids, from = 'cid', to = 'property', properties= properties
    )]


treatment_annotations <- merge(
    compound_nameToCIDS, 
    pubchemProperties[, cids := as.character(CID)], 
    by.x= "cids",  by.y = "cids", all.x = TRUE)


saveRDS(treatment_annotations, file = OUTPUT$treatmentMetadata)