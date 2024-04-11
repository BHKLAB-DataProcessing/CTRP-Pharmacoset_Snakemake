## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
# This snippet is run at the beginning of a snakemake run to setup the env
# Helps to load the workspace if the script is run independently or debugging
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

    # Assuming that this script is named after the rule
    # Saves the workspace to "resources/"annotate_sampleMetadata"
    file.path("resources", paste0(snakemake@rule, ".RData")) |> 
        save.image()
}else{
    # If the snakemake object does not exist, load the workspace
    file.path("resources", "annotate_sampleMetadata.RData") |>
        load()
}


###############################################################################
# Load INPUT
###############################################################################
sampleMetadata <- data.table::fread(INPUT$rawSampleMetadata, sep = "\t", header = TRUE)
str(sampleMetadata)
message("Number of samples: ", nrow(sampleMetadata))

# Set this for mapping 
options("mc.cores" = THREADS)

options("log_level" = "INFO")   # AnnotationGx logging level

###############################################################################
# Main Script
###############################################################################
# Using the depMapID over the cell line name since this ensures the exact match
# that depmap would use
message("Mapping cell line names to Cellosaurus accessions")
options("log_level" = "INFO")
(mapped_cells <- AnnotationGx::mapCell2Accession(
    ids = sampleMetadata[sampleid != "", sampleid]
))
mapped_cells <- unique(mapped_cells[!is.na(accession),])
mapped_cells <- merge(sampleMetadata, mapped_cells, by.x = "sampleid", by.y = "query", all.x = FALSE)
message("Successful mappings: ", nrow(mapped_cells))


failed <- sampleMetadata[!sampleid %in% mapped_cells$sampleid, ]

message("Number of failed mappings: ", nrow(failed))
print(missing <- failed)




# Here we parse the sample id on the _ character and use the first part to map to the cellosaurus
# accession number. Also use fuzzy matching to try and get the correct cell line name
message("Trying again using the sample name")

missing[, c("cellLineName", "accession") := {
    result <- AnnotationGx::mapCell2Accession(sampleid, fuzzy = T)
    list(result$cellLineName, result$accession)
}]
missing


# Combine the mapped and missing data
(annotated_sampleMetadata <- data.table::rbindlist(
    list(
        mapped_cells, 
        missing
    ), fill = TRUE
)[order(sampleid)])

# Rename columns from ctrp
data.table::setnames(
    annotated_sampleMetadata,
    c("tissueid", "master_ccl_id", "sample_availability", "ccle_primary_hist"),
    c("CTRP.tissueid", "CTRP.master_ccl_id", "CTRP.sample_availability", "CTRP.ccle_primary_hist"),
    skip_absent = TRUE
)

# Rename the columns from cellosaurus
data.table::setnames(
    annotated_sampleMetadata,
    c("cellLineName", "accession"),
    c("cellosaurus.cellLineName", "cellosaurus.accession")
)

message("Annotating sample metadata")

annotated_accessions <- AnnotationGx::annotateCellAccession(
    accessions = annotated_sampleMetadata[!is.na(cellosaurus.accession), cellosaurus.accession],
)
message("Number of annotated accessions: ", nrow(annotated_accessions))

message("Number of unique categories: ")
annotated_accessions[, .N, by = "category"]

message("Number of unique sexOfCell: ")
annotated_accessions[, .N, by = "sexOfCell"]


annotated_accessions[, synonyms := sapply(synonyms, function(x) paste(x, collapse = "; "))]
annotated_accessions[, diseases := sapply(diseases, function(x) paste(x, collapse = "; "))]


annotated_accessions[, c("crossReferences", "hierarchy", "comments") := NULL]
names(annotated_accessions) <- paste0("cellosaurus.", names(annotated_accessions))

annotated_accessions <- unique(annotated_accessions)



## ------------------------------------------------------------------------- ##


# merge
final_annotated <- merge(
    annotated_sampleMetadata, 
    annotated_accessions,
    by= c("cellosaurus.accession", "cellosaurus.cellLineName"),
    all.x = TRUE
) |> unique()

data.table::setcolorder(
    final_annotated,
    c("cellosaurus.accession", "cellosaurus.cellLineName"),
    before= "cellosaurus.category" 
)

final_annotated <- final_annotated[order(sampleid)]
str(final_annotated)



###############################################################################
# Save OUTPUT 
###############################################################################
message("Saving sampleMetadata to: ", OUTPUT$sampleMetadata)
dir.create(dirname(OUTPUT$sampleMetadata), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(
    final_annotated, 
    file = OUTPUT$sampleMetadata, 
    quote = TRUE, 
    sep = "\t", 
    na = "NA", 
    col.names = TRUE
)
