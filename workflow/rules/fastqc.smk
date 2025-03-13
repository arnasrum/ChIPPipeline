rule fastqc:
    input:
        raw = RESOURCES + "/reads/{sample}.fastq.gz",
        trimmed = RESULTS + "/{tool}/{sample}.fastq.gz",
    output:
        raw = multiext(RESULTS + "/fastqc/{tool}/raw/{sample}_fastqc.", "zip", "html"),
        trimmed = multiext(RESULTS + "/fastqc/{tool}/trimmed/{sample}_fastqc.","zip","html")
    params:
        outputPath = lambda wildcards: f"{RESULTS}/fastqc/{wildcards.tool}"
    conda:
        "../envs/fastqc.yml"
    log:
        config['logs_path'] + "/fastqc/{tool}/{sample}.log",
    threads:
        int(config["fastqc"]["threads"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb = int(config["fastqc"]["mem_mb"])
    shell:
        """
        exec > {log} 2>&1
        fastqc --dir {resources.tmpdir} -t {threads} -o {params.outputPath}/raw --memory {resources.mem_mb} {input.raw} 
        fastqc --dir {resources.tmpdir} -t {threads} -o {params.outputPath}/trimmed --memory {resources.mem_mb} {input} 
        """



