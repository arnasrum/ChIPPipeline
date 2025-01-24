
# Add prefix to these paths
rule fastqc_after_trimming:
    input:
        config['results_path'] + "/" + "{tool}/{sample}.fastq",
    output:
        multiext(config['results_path'] + "/fastqc/{tool}/{sample}_fastqc.", "zip", "html")
    params:
        outputPath = lambda wildcards: f"{config['results_path']}/fastqc/" + wildcards.tool
    conda:
        "../envs/fastqc.yml"
    log:
        config['logs_path'] + "/fastqc/{tool}/{sample}.log"
    threads:
        int(config["fastqc"]["threads"])
    resources:
        tmpdir=config["temp_path"],
        mem_mb=1024
    shell:
        """
        exec > {log} 2>&1
        fastqc -t {threads} -o {params.outputPath} --memory {resources.mem_mb} {input} 
        """

rule fastqc_before_trimming:
    input:
        RESOURCES + "/reads/{sample}.fastq",
    output:
        multiext(config['results_path'] + "/fastqc/raw/{sample}_fastqc.", "zip", "html")
    params:
        outputPath = lambda wildcards: f"{config['results_path']}/fastqc/raw"
    conda:
        "../envs/fastqc.yml"
    log:
        config['logs_path'] + "/fastqc/raw/{sample}.log"
    threads:
        int(config["fastqc"]["threads"])
    resources:
        tmpdir=config["temp_path"],
        mem_mb=1024
    shell:
        """
        exec > {log} 2>&1
        fastqc -t {threads} -o {params.outputPath} --memory {resources.mem_mb} {input} 
        """



