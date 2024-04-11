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
    file.path("resources", "annotate_treatmentMetadata.RData") |>
        load()
}

#############################################################################
# Load INPUT
#############################################################################

# zipDir <- file.path(RESOURCES$tmpdir, snakemake@rule, tools::file_path_sans_ext(INPUT$tr))
# dir.create(zipDir, recursive = TRUE, showWarnings = FALSE)
# print(paste0("Unzipping ", INPUT$tr, " into ", zipDir))
# unzip(INPUT$tr, exdir = zipDir)
# files <- paste0(zipDir, "/", list.files(zipDir))


#############################################################################
# Main Script
#############################################################################
message("Reading in treatment metadata")
treatmentMetadata <- data.table::fread(INPUT$rawTreatmentMetadata)

ann_cols <- treatmentMetadata[, .(cpd_name, broad_cpd_id)]

message("Mapping compound names to PubChem CIDs")
(compound_nameToCIDS <- AnnotationGx::mapCompound2CID(ann_cols[, cpd_name], first = TRUE))

properties=c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES')
message(
    "Getting the following properties from PubChem: ", 
    paste(properties, collapse= " "), " for ", nrow(compound_nameToCIDS), " compounds")

pubchemProperties <- AnnotationGx::getPubchemCompound(
        ids = compound_nameToCIDS[!is.na(cids), cids], 
        from = 'cid', 
        to = 'property', 
        properties = properties
    )


treatment_annotations <- merge(
    compound_nameToCIDS[, cids := as.character(cids)], 
    pubchemProperties[, cids := as.character(CID)], 
    by.x= "cids",  by.y = "cids", all.x = TRUE)

names(treatment_annotations) <- paste0("pubchem.", names(treatment_annotations))

show(treatment_annotations)

# Unichem Mapping
# _______________
message("Getting Unichem sources")
unichem_sources <- AnnotationGx::getUnichemSources(T)
str(unichem_sources)

sources_of_interest <- c("chembl", "drugbank", "chebi", "phamgkb", "lincs", "clinicaltrials", "nih_ncc", "fdasrs", "pharmgkb", "rxnorm")
sourceID <- unichem_sources[Name == "pubchem", SourceID]

message("Mapping PubChem CIDs to other databases")

result <- parallel::mclapply(
    treatment_annotations[!is.na(pubchem.cids), pubchem.cids], 
    function(cid){
        tryCatch({
            # query the Unichem database for the compound
            result <- AnnotationGx::queryUnichemCompound(type = "sourceID", compound = cid, sourceID = sourceID)

            # subset the results to only the sources of interest
            subset <- result$External_Mappings[Name %in% sources_of_interest, .(compoundID, Name)]
            
            # make `Name` = the column names and the `values` =  compoundID 
            subset$cid <- cid # need to add column using `cid` so we can merge back later
            data.table::dcast(subset, cid ~ Name, value.var = "compoundID", fun.aggregate = list)
        }, error = function(e) NULL)
    }, mc.cores = 30
  ) |> data.table::rbindlist(fill = T)

unichem_mappings <- data.table::copy(result)
# for each column, if its a list then make it a string with a comma separator
for(col in names(unichem_mappings)){
  if(is.list(unichem_mappings[[col]])){
    unichem_mappings[[col]] <- sapply(unichem_mappings[[col]], function(x) paste(x, collapse = ","))
  }
}

data.table::setkey(unichem_sources, Name)

# Get the corresponding source name for each column and tag with unichem. (i.e chembl -> unichem.ChEMBL)
# Note `cid` becomes unichem.NA because `cid` isnt in unichem_sources, which is fine since we lose that column
names(unichem_mappings) <- paste("unichem", unichem_sources[names(unichem_mappings), gsub(" ", "_", NameLabel)], sep = ".")

show(unichem_mappings)

unichem_annotated_treatmentMetadata <- merge(
    treatment_annotations, 
    unichem_mappings, 
    by.x = "pubchem.cids", by.y = "unichem.NA", all.x = TRUE
)

message("Using unichem.chembl to map to chembl Database")

has_chembl_ids <- unichem_annotated_treatmentMetadata[!is.na(unichem.ChEMBL) & unichem.ChEMBL != "", ]

chembl_mechanisms <- parallel::mclapply(
    has_chembl_ids$unichem.ChEMBL, 
    AnnotationGx::getChemblMechanism,
    mc.cores = 30
) 

chembl_mechanisms <- data.table::rbindlist(chembl_mechanisms, fill = T)

dt <- chembl_mechanisms[!is.na(molecule_chembl_id),]

# remove all columns in dt that are empty
dt <- dt[, .SD, .SDcols = colSums(is.na(dt)) < nrow(dt)]

chembl_dt <- dt[
    , 
    .(
        action_type, disease_efficacy, max_phase, mec_id, mechanism_of_action, 
        molecule_chembl_id, target_chembl_id, site_id, parent_molecule_chembl_id
    )
]

names(chembl_dt) <- paste("chembl", names(chembl_dt), sep = ".")


all_annotated_treatmentMetadata_chembl <- merge(
  unique(unichem_annotated_treatmentMetadata), 
  unique(chembl_dt),
  by.x = "unichem.ChEMBL", 
  by.y = "chembl.molecule_chembl_id", 
  all.x = T, 
  incomparables = NA_character_
)

# Renaming the `treatmentMetadata` columns before merging
# _______________________________________________________
names(treatmentMetadata) <- paste("CTRP.", names(treatmentMetadata), sep = "")
# rename "CTRP.treatmentid" to just "treatmentid"
data.table::setnames(treatmentMetadata, "CTRP.treatmentid", "treatmentid")
str(treatmentMetadata)

final_treatmentMetadata <- 
    merge(
        treatmentMetadata, 
        all_annotated_treatmentMetadata_chembl, 
        by.x = "treatmentid", 
        by.y = "pubchem.name", 
        all.x = TRUE
    )

# reorder some of the columns here for ease:

data.table::setcolorder(final_treatmentMetadata, "treatmentid")

data.table::setcolorder(final_treatmentMetadata, "unichem.ChEMBL", before = "chembl.action_type")

names(final_treatmentMetadata)

## --------------------- Save OUTPUT ------------------- ##
dir.create(dirname(OUTPUT$treatmentMetadata), showWarnings = FALSE, recursive = TRUE)

data.table::fwrite(
    final_treatmentMetadata, 
    OUTPUT$treatmentMetadata,
    sep = "\t",
    quote = FALSE)

