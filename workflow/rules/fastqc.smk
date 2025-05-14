rule fastqc_unprocessed:
    input:
        raw = f"{RESOURCES}/reads/{{sample}}.fastq.gz",
    output:
        raw = multiext(f"{RESULTS}/fastqc/unprocessed/{{sample}}_fastqc.", "zip", "html"),
    params:
        outputPath = lambda wildcards: f"{RESULTS}/fastqc/unprocessed"
    conda:
        "../envs/fastqc.yml"
    log:
        f"{LOGS}/fastqc/{{sample}}.log",
    threads:
        int(config["fastqc"]["threads"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: int(config["fastqc"]["mem_mb"]) * attempt,
        runtime = lambda wildcards, attempt: int(config["fastqc"]["runtime"]) * attempt,
    shell:
        """
        exec > {log} 2>&1
        fastqc --dir {resources.tmpdir} -t {threads} -o {params.outputPath} {input.raw} 
        """

rule fastqc_trimmed:
    input:
        trimmed = f"{RESULTS}/{{tool}}/{{sample}}.fastq.gz",
    output:
        trimmed = multiext(f"{RESULTS}/fastqc/{{tool}}/{{sample}}_fastqc.","zip","html")
    params:
        outputPath = lambda wildcards: f"{RESULTS}/fastqc/{wildcards.tool}"
    conda:
        "../envs/fastqc.yml"
    log:
        f"{LOGS}/fastqc/{{tool}}/{{sample}}.log",
    threads:
        int(config["fastqc"]["threads"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: int(config["fastqc"]["mem_mb"]) * attempt,
        runtime = lambda wildcards, attempt: int(config["fastqc"]["runtime"]) * attempt,
    shell:
        """
        exec > {log} 2>&1
        fastqc --dir {resources.tmpdir} -t {threads} -o {params.outputPath} {input.trimmed} 
        """
