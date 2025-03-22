rule trim_galore:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/trim_galore/{sample}{extension}", extension=fastq_file_extensions, allow_missing=True))
    conda:
        "../envs/trim.yml"
    params:
        args = config["trim_galore"]["args"],
        output_dir = RESULTS + "/trim_galore",
    threads:
        int(config["trim_galore"]["threads"])
    log:
        LOGS + "/trim_galore/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/trim_galore/{sample}.txt", config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['trim_galore']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['trim_galore']['runtime'] * attempt
    script:
        "../scripts/tools/trim_galore.py"

rule cutadapt:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/cutadapt/{sample}{extension}", extension=fastq_file_extensions, allow_missing=True))
    conda:
        "../envs/cutadapt.yml"
    params:
        args = config["cutadapt"]["args"],
    threads:
        int(config["cutadapt"]["threads"])
    log:
        LOGS + "/cutadapt/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/cutadapt/{sample}.txt", config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['cutadapt']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['cutadapt']['runtime'] * attempt
    script:
        "../scripts/tools/cutadapt.py"

rule fastp:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/fastp/{sample}{extension}", extension=fastq_file_extensions, allow_missing=True))
    conda:
        "../envs/trim.yml"
    params:
        args = config["fastp"]["args"],
    threads:
        int(config["fastp"]["threads"])
    log:
        LOGS + "/fastp/{sample}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb= lambda wildcards,attempt: config['fastp']['mem_mb'] * attempt,
        runtime= lambda wildcards,attempt: config['fastp']['runtime'] * attempt
    benchmark:
        repeat(BENCHMARKS + "/fastp/{sample}.txt", config["benchmark_repeat_trim"])
    script:
        "../scripts/tools/fastp.py"

rule trimmomatic:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/trimmomatic/{sample}{extension}", extension=fastq_file_extensions, allow_missing=True))
    conda:
        "../envs/trim.yml"
    params:
        args = config["trimmomatic"]["args"],
        run_options = config['trimmomatic']['run_options']
    threads:
        int(config["trimmomatic"]["threads"])
    log:
        LOGS + "/trimmomatic/{sample}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=config['trimmomatic']['mem_mb'],
        runtime=config['trimmomatic']['runtime']
    benchmark:
        repeat(BENCHMARKS + "/trimmomatic/{sample}.txt", config["benchmark_repeat_trim"])
    script:
        "../scripts/tools/trimmomatic.py"

