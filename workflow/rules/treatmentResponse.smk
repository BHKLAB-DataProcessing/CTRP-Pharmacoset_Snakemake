configfile: "config/pipeline.yaml"


rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("workflow/scripts")
# create a Path object to standardize paths to treatmentResponse data
trDir = Path("treatmentResponse")


rule build_treatmentResponseExperiment:
    input:
        tr=rawdata / trDir / "CTRPv{release}.zip",
        treatmentMetadata=procdata / metadata / "CTRPv{release}_treatmentMetadata.tsv",
        sampleMetadata=procdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
    output:
        tre=procdata / trDir / "CTRPv{release}_treatmentResponseExperiment.RDS",
        sensInfo=procdata / trDir / "CTRPv{release}_sensitivityInfo.RDS",
    log:
        logs / "{release}" / "build_treatmentResponseExperiment.log",
    threads: 30
    script:
        ".." / scripts / trDir / "build_treatmentResponseExperiment.R"
