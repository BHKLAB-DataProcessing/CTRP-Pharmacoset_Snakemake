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

rule all:
    input:
        expand(
            procdata
            / "treatmentResponse"
            / "CTRPv{release}_treatmentResponseExperiment.RDS",
            release=["2.0"],
        ),
