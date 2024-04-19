from pathlib import Path

configfile: "config/pipeline.yaml"

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("workflow/scripts")

include: "workflow/rules/download.smk"
include: "workflow/rules/metadata.smk"
include: "workflow/rules/treatmentResponse.smk"

rule build_PharmacoSet:
    input:
        treatmentResponseExperiment = expand(
            procdata / "treatmentResponse" / "CTRPv{release}_treatmentResponseExperiment.RDS",
            release=["2.0"],
        ),
        multiAssayExperiment = expand(
            procdata / "CTRPv{release}_multiAssayExperiment.RDS",
            release=["2.0"],
        ),
        treatmentMetadata = expand(
            procdata / metadata / "CTRPv{release}_treatmentMetadata.tsv",
            release=["2.0"],
        ),
        sampleMetadata = expand(
            procdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
            release=["2.0"],
        ),

rule build_MultiAssayExperiment:
    input:
        sampleMetadata = 
            procdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
    output:
        procdata / "CTRPv{release}_multiAssayExperiment.RDS",
    log:
        logs / "build_MultiAssayExperiment_{release}.log",
    conda:
        "workflow/envs/PharmacoSet.yaml"
    threads:
        1,
    script:
        scripts / "build_MultiAssayExperiment.R"