

rule fastqc_after_trim:
    input:
        f"results/{config["trimmer"]}/{{id}}.fastq",
    output:
        expand("results/fastqc/{trimmer}/{id}_fastqc.{ext}", trimmer=config["trimmer"], ext=["zip", "html"], allow_missing=True)
    params:
        outputPath = "results/fastqc/" + config["trimmer"]
    shell:
        """
        fastqc -o {params.outputPath} {input} 
        """


