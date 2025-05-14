ruleorder: trim_galore_PE > trim_galore_SE
ruleorder: fastp_PE > fastp_SE
ruleorder: cutadapt_PE > cutadapt_SE
ruleorder: trimmomatic_PE > trimmomatic_SE

rule trim_galore_SE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        temp(f"{RESULTS}/trim_galore/{{sample}}.fastq.gz")
    conda:
        "../envs/trim.yml"
    params:
        args = config["trim_galore"]["args"],
        output_dir = f"{RESULTS}/trim_galore",
    threads:
        int(config["trim_galore"]["threads"])
    log:
        f"{LOGS}/trim_galore/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/trim_galore/{{sample}}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['trim_galore']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['trim_galore']['runtime'] * attempt
    script:
        "../scripts/tools/trim_galore.py"

rule trim_galore_PE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        temp(f"{RESULTS}/trim_galore/{{sample}}_1.fastq.gz"),
        temp(f"{RESULTS}/trim_galore/{{sample}}_2.fastq.gz"),
    conda:
        "../envs/trim.yml"
    params:
        args=config["trim_galore"]["args"],
        output_dir=f"{RESULTS}/trim_galore",
    threads:
        int(config["trim_galore"]["threads"])
    log:
        f"{LOGS}/trim_galore/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/trim_galore/{{sample}}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['trim_galore']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['trim_galore']['runtime'] * attempt
    script:
        "../scripts/tools/trim_galore.py"

rule cutadapt_SE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        temp(f"{RESULTS}/cutadapt/{{sample}}.fastq.gz")
    conda:
        "../envs/trim.yml"
    params:
        args = config["cutadapt"]["args"],
    threads:
        int(config["cutadapt"]["threads"])
    log:
        f"{LOGS}/cutadapt/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/cutadapt/{{sample}}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['cutadapt']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['cutadapt']['runtime'] * attempt
    script:
        "../scripts/tools/cutadapt.py"

rule cutadapt_PE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        temp(f"{RESULTS}/cutadapt/{{sample}}_1.fastq.gz"),
        temp(f"{RESULTS}/cutadapt/{{sample}}_2.fastq.gz"),
    conda:
        "../envs/trim.yml"
    params:
        args=config["cutadapt"]["args"],
    threads:
        int(config["cutadapt"]["threads"])
    log:
        f"{LOGS}/cutadapt/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/cutadapt/{{sample}}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['cutadapt']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['cutadapt']['runtime'] * attempt
    script:
        "../scripts/tools/cutadapt.py"

rule fastp_SE:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        temp(f"{RESULTS}/fastp/{{sample}}.fastq.gz")
    conda:
        "../envs/trim.yml"
    params:
        args = config["fastp"]["args"],
    threads:
        int(config["fastp"]["threads"])
    log:
        f"{LOGS}/fastp/{{sample}}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb= lambda wildcards,attempt: config['fastp']['mem_mb'] * attempt,
        runtime= lambda wildcards,attempt: config['fastp']['runtime'] * attempt
    benchmark:
        f"{BENCHMARKS}/fastp/{{sample}}.txt"
    script:
        "../scripts/tools/fastp.py"

rule fastp_PE:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        temp(f"{RESULTS}/fastp/{{sample}}_1.fastq.gz"),
        temp(f"{RESULTS}/fastp/{{sample}}_2.fastq.gz")
    conda:
        "../envs/trim.yml"
    params:
        args=config["fastp"]["args"],
    threads:
        int(config["fastp"]["threads"])
    log:
        f"{LOGS}/fastp/{{sample}}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['fastp']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['fastp']['runtime'] * attempt
    benchmark:
        f"{BENCHMARKS}/fastp/{{sample}}.txt"
    script:
        "../scripts/tools/fastp.py"


rule trimmomatic_SE:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        temp(f"{RESULTS}/trimmomatic/{{sample}}.fastq.gz")
    conda:
        "../envs/trim.yml"
    params:
        args = config["trimmomatic"]["args"],
        trimming_steps = config['trimmomatic']['trimming_steps']
    threads:
        int(config["trimmomatic"]["threads"])
    log:
        f"{LOGS}/trimmomatic/{{sample}}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=config['trimmomatic']['mem_mb'],
        runtime=config['trimmomatic']['runtime']
    benchmark:
        f"{BENCHMARKS}/trimmomatic/{{sample}}.txt"
    script:
        "../scripts/tools/trimmomatic.py"

rule trimmomatic_PE:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        temp(f"{RESULTS}/trimmomatic/{{sample}}_1.fastq.gz"),
        temp(f"{RESULTS}/trimmomatic/{{sample}}_2.fastq.gz")
    conda:
        "../envs/trim.yml"
    params:
        args = config["trimmomatic"]["args"],
        trimming_steps = config['trimmomatic']['trimming_steps']
    threads:
        int(config["trimmomatic"]["threads"])
    log:
        f"{LOGS}/trimmomatic/{{sample}}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=config['trimmomatic']['mem_mb'],
        runtime=config['trimmomatic']['runtime']
    benchmark:
        f"{BENCHMARKS}/trimmomatic/{{sample}}.txt"
    script:
        "../scripts/tools/trimmomatic.py"
