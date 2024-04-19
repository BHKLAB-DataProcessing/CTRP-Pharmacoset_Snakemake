from pathlib import Path

metadata = Path(config["directories"]["metadata"])

annotationGx_Docker = config["containers"]["annotationGx"]

rule preprocess_Metadata:
    input:
        tr=rawdata / "treatmentResponse" / "CTRPv{release}.zip",
    output:
        rawSampleMetadata=rawdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
        rawTreatmentMetadata=rawdata / metadata / "CTRPv{release}_treatmentMetadata.tsv",
    log:
        logs / "CTRPv{release}_preprocess_Metadata.log",
    container:
        annotationGx_Docker,
    script:
        "../scripts" / metadata / "preprocess_Metadata.R"


rule annotate_treatmentMetadata:
    input:
        rawTreatmentMetadata=rawdata / metadata / "CTRPv{release}_treatmentMetadata.tsv",
    output:
        treatmentMetadata=procdata / metadata / "CTRPv{release}_treatmentMetadata.tsv",
    log:
        logs / "CTRPv{release}_annotate_treatmentMetadata.log",
    container:
        annotationGx_Docker,
    threads:
        4
    script:
       "../scripts" /  metadata / "annotate_treatmentMetadata.R"


rule annotate_sampleMetadata:
    input:
        rawSampleMetadata=rawdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
    output:
        sampleMetadata=procdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
    log:
        logs / "CTRPv{release}_annotate_sampleMetadata.log",
    container:
        annotationGx_Docker,
    threads:
        4
    script:
        "../scripts"/ metadata / "annotate_sampleMetadata.R"
