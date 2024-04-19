from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

treatmentResponse = config["treatmentResponse"]

rule download_treatmentResponse_AND_metadata:
    input:
        lambda wc: HTTP.remote(
            treatmentResponse[wc.release]["url"],
        )
    output:
        tr = rawdata / "treatmentResponse" / "CTRPv{release}.zip" 
    shell:
        """
        mv {input[0]} {output.tr}
        """