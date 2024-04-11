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




#############################################################################
# Main Script




## ----------------------------------------------------- ##




## --------------------- Save OUTPUT ------------------- ##