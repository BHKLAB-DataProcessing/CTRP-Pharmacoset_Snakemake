from pathlib import Path


configfile: "workflow/config/pipeline.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("workflow/scripts")

treatmentResponse = config["treatmentResponse"]

storage:
    provider = "http"


rule all:
    input:
        expand(
            procdata / "treatmentResponse" / "CTRPv{release}_treatmentResponseExperiment.RDS",
                release=["2.0"]
        ),
        expand(
            procdata / metadata / "CTRPv{release}_treatmentMetadata.RDS"
                release=["2.0"]
        ),

rule download_treatmentResponse:
    input:
        lambda wc: storage.http(
            treatmentResponse[wc.release]["url"],
        )
    output:
        tr = rawdata / "treatmentResponse" / "CTRPv{release}.zip" 
    shell:
        # echo $(basename {input[0]} .zip)
        """
        mv {input[0]} {output.tr}
        """

rule build_treatmentResponseExperiment:
    input:
        tr = rawdata / "treatmentResponse" / "CTRPv{release}.zip"
    output:
        tre = procdata / "treatmentResponse" / "CTRPv{release}_treatmentResponseExperiment.RDS",
        sensInfo = procdata / "treatmentResponse" / "CTRPv{release}_sensitivityInfo.RDS"
    log:
        logs / "{release}" / "build_treatmentResponseExperiment.log"
    threads: 
        30
    script:
        scripts / "treatmentResponse" / "build_treatmentResponseExperiment.R"

rule annotate_treatmentMetadata:
    input:
        tr = rawdata / "treatmentResponse" / "CTRPv{release}.zip",
    output:
        treatmentMetadata = procdata / metadata / "CTRPv{release}_treatmentMetadata.RDS"
    script:
        scripts / metadata / "preprocess_treatmentMetadata.R"

