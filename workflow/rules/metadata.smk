
rule preprocess_Metadata:
    input:
        tr=rawdata / "treatmentResponse" / "CTRPv{release}.zip",
    output:
        rawSampleMetadata=rawdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
        rawTreatmentMetadata=rawdata / metadata / "CTRPv{release}_treatmentMetadata.tsv",
    log:
        logs / "CTRPv{release}_preprocess_Metadata.log",
    script:
        ".." / scripts / metadata / "preprocess_Metadata.R"


rule annotate_treatmentMetadata:
    input:
        rawTreatmentMetadata=rawdata / metadata / "CTRPv{release}_treatmentMetadata.tsv",
    output:
        treatmentMetadata=procdata / metadata / "CTRPv{release}_treatmentMetadata.tsv",
    script:
        ".." / scripts / metadata / "annotate_treatmentMetadata.R"


rule annotate_sampleMetadata:
    input:
        rawSampleMetadata=rawdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
    output:
        sampleMetadata=procdata / metadata / "CTRPv{release}_sampleMetadata.tsv",
    script:
        ".." / scripts / metadata / "annotate_sampleMetadata.R"
