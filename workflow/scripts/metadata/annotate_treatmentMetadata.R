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
    file.path("resources", "preprocess_treatmentMetadata.RData") |>
        load()
}

#############################################################################
# Load INPUT
#############################################################################

zipDir <- file.path(RESOURCES$tmpdir, snakemake@rule, tools::file_path_sans_ext(INPUT$tr))
dir.create(zipDir, recursive = TRUE, showWarnings = FALSE)
print(paste0("Unzipping ", INPUT$tr, " into ", zipDir))
unzip(INPUT$tr, exdir = zipDir)
files <- paste0(zipDir, "/", list.files(zipDir))


#############################################################################
# Main Script



# Get the treatment metadata
# --------------------------
treatmentMetadata <- 
    files[grepl(pattern = "*per_compound*", x = files)] |>
    data.table::fread() 
treatmentMetadata[, treatmentid := cpd_name]

ann_cols <- treatmentMetadata[, .(cpd_name, broad_cpd_id)]


(compound_nameToCIDS <- AnnotationGx::mapCompound2CID(ann_cols[, cpd_name], first = TRUE))

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



# annotations <- c('ChEMBL ID', 'NSC Number', 'Drug Induced Liver Injury', 'CAS', 'ATC Code')

# message("Annotating with ChEMBL ID")
# treatment_annotations[, 'pubchem.ChEMBL.ID' := lapply(cids, function(x) AnnotationGx::annotatePubchemCompound(cid = x, heading = 'ChEMBL ID')), by = cids]

# message("Annotating with ChEMBL ID")
# pubchem_annotated[, 'pubchem.ChEMBL.ID' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'ChEMBL ID')), by = pubchem.CID]

# message("Annotating with NSC Number")
# pubchem_annotated[, 'pubchem.NSC.Number' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'NSC Number')), by = pubchem.CID]

# message("Annotating with Drug Induced Liver Injury")
# pubchem_annotated[, 'pubchem.DILI.Status' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'Drug Induced Liver Injury')), by = pubchem.CID]

# message("Annotating with CAS Number")
# pubchem_annotated[, 'pubchem.CAS.Number' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'CAS')), by = pubchem.CID]

# message("Annotating with ATC Code")
# pubchem_annotated[, 'pubchem.ATC.Code' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'ATC Code')), by = pubchem.CID]

# message("Annotating with Synonyms")
# annotated_treatments <- merge(treatmentMetadata, pubchem_annotated, by.x = "GDSC.treatmentid", by.y = "GDSC.treatmentid", all.x = TRUE)




## ----------------------------------------------------- ##




## --------------------- Save OUTPUT ------------------- ##

data.table::fwrite(
    treatment_annotations, 
    OUTPUT$treatmentMetadata,
    sep = "\t",
    quote = FALSE)