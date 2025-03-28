ruleorder: trim_galore_SE > trim_galore_PE
ruleorder: fastp_SE > fastp_PE
ruleorder: cutadapt_SE > cutadapt_PE
ruleorder: trimmomatic_SE > trimmomatic_PE

rule trim_galore_SE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(RESULTS + "/trim_galore/{sample}.fastq.gz")
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

rule trim_galore_PE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(RESULTS + "/trim_galore/{sample}_1.fastq.gz"),
        temp(RESULTS + "/trim_galore/{sample}_2.fastq.gz"),
    conda:
        "../envs/trim.yml"
    params:
        args=config["trim_galore"]["args"],
        output_dir=RESULTS + "/trim_galore",
    threads:
        int(config["trim_galore"]["threads"])
    log:
        LOGS + "/trim_galore/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/trim_galore/{sample}.txt",config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['trim_galore']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['trim_galore']['runtime'] * attempt
    script:
        "../scripts/tools/trim_galore.py"


rule cutadapt_SE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(RESULTS + "/cutadapt/{sample}.fastq.gz")
    conda:
        "../envs/trim.yml"
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

rule cutadapt_PE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(RESULTS + "/cutadapt/{sample}_1.fastq.gz"),
        temp(RESULTS + "/cutadapt/{sample}_2.fastq.gz"),
    conda:
        "../envs/cutadapt.yml"
    params:
        args=config["cutadapt"]["args"],
    threads:
        int(config["cutadapt"]["threads"])
    log:
        LOGS + "/cutadapt/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/cutadapt/{sample}.txt",config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['cutadapt']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['cutadapt']['runtime'] * attempt
    script:
        "../scripts/tools/cutadapt.py"

rule fastp_SE:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(RESULTS + "/fastp/{sample}.fastq.gz")
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

rule fastp_PE:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(RESULTS + "/fastp/{sample}_1.fastq.gz"),
        temp(RESULTS + "/fastp/{sample}_2.fastq.gz")
    conda:
        "../envs/trim.yml"
    params:
        args=config["fastp"]["args"],
    threads:
        int(config["fastp"]["threads"])
    log:
        LOGS + "/fastp/{sample}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['fastp']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['fastp']['runtime'] * attempt
    benchmark:
        repeat(BENCHMARKS + "/fastp/{sample}.txt",config["benchmark_repeat_trim"])
    script:
        "../scripts/tools/fastp.py"


rule trimmomatic_SE:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(RESULTS + "/trimmomatic/{sample}.fastq.gz"),
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

rule trimmomatic_PE:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(RESULTS + "/trimmomatic/{sample}_1.fastq.gz"),
        temp(RESULTS + "/trimmomatic/{sample}_2.fastq.gz"),
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
