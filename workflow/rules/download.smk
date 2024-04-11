
treatmentResponse = config["treatmentResponse"]

storage:
    provider = "http"

rule download_treatmentResponse_AND_metadata:
    input:
        lambda wc: storage.http(
            treatmentResponse[wc.release]["url"],
        )
    output:
        tr = rawdata / "treatmentResponse" / "CTRPv{release}.zip" 
    shell:
        """
        mv {input[0]} {output.tr}
        """