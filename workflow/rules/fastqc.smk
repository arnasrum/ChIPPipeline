

rule fastqc_after_trim:
    input:
        f"results/{config['trimmer']}/{{sample}}.fastq",
    output:
        multiext(f"results/fastqc/{config['trimmer']}/{{sample}}.", "zip", "html")
    params:
        outputPath = "results/fastqc/" + config["trimmer"]
    shell:
        """
        fastqc -o {params.outputPath} {input} 
        """


