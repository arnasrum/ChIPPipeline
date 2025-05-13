ruleorder: get_fastq_PE > get_fastq_SE
ruleorder: concatenate_runs_PE > concatenate_runs_SE

rule get_fastq_PE:
    output:
        temp(f"{RESOURCES}/reads/{{accession}}_1.fastq"),
        temp(f"{RESOURCES}/reads/{{accession}}_2.fastq")
    log:
        f"{LOGS}/fastq-dump_pe/{{accession}}.log"
    params:
        extra=str(config["fastq-dump"]["args"])
    threads:
        int(config["fastq-dump"]["threads"])
    resources:
        tmpdir = TEMP,
        cpus_per_task = int(config["fastq-dump"]["threads"]),
        runtime = int(config["fastq-dump"]["runtime"]),
        mem_mb = int(config["fastq-dump"]["mem_mb"])
    wrapper:
        "v5.10.0/bio/sra-tools/fasterq-dump"

rule get_fastq_SE:
    output:
        temp(f"{RESOURCES}/reads/{{accession}}.fastq")
    log:
        f"{LOGS}/fastq-dump_se/{{accession}}.log"
    params:
        extra=lambda wildcards, resources: str(config["fastq-dump"]["args"])
    threads:
        int(config["fastq-dump"]["threads"])
    resources:
        tmpdir = TEMP,
        cpus_per_task = int(config["fastq-dump"]["threads"]),
        runtime = int(config["fastq-dump"]["runtime"]),
        mem_mb = int(config["fastq-dump"]["mem_mb"])
    wrapper:
        "v5.10.0/bio/sra-tools/fasterq-dump"

rule concatenate_runs_SE:
    input:
        lambda wildcards: concatenate_runs_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        f"{RESOURCES}/reads/{{sample}}.fastq.gz"
    log:
        f"{LOGS}/concatenate/{{sample}}.log"
    resources:
        tmpdir=TEMP
    conda:
        "../envs/input.yml"
    script:
        "../scripts/concatenate_runs.py"

rule concatenate_runs_PE:
    input:
        lambda wildcards: concatenate_runs_input(wildcards.sample, RESOURCES, pipeline_config)
    output:
        f"{RESOURCES}/reads/{{sample}}_1.fastq.gz",
        f"{RESOURCES}/reads/{{sample}}_2.fastq.gz"
    log:
        f"{LOGS}/concatenate/{{sample}}.log"
    resources:
        tmpdir=TEMP
    conda:
        "../envs/input.yml"
    script:
        "../scripts/concatenate_runs.py"

rule handle_provided_samples_SE:
    input:
        lambda wildcards: symlink_input(wildcards.file_name, pipeline_config)["read1"]["path"],
    output:
        f"{RESOURCES}/reads/{{file_name}}.fastq.gz"
    log:
        f"{LOGS}/provided_files/{{file_name}}.log"
    resources:
        tmpdir=TEMP
    conda:
        "../envs/input.yml"
    script:
        "../scripts/handle_provided_files.py"

rule handle_provided_samples_PE:
    input:
        lambda wildcards: symlink_input(wildcards.file_name, pipeline_config)["read1"]["path"],
        lambda wildcards: symlink_input(wildcards.file_name, pipeline_config)["read2"]["path"]
    output:
        f"{RESOURCES}/reads/{{file_name}}_1.fastq.gz",
        f"{RESOURCES}/reads/{{file_name}}_2.fastq.gz"
    log:
        f"{LOGS}/provided_files/{{file_name}}.log"
    resources:
        tmpdir=TEMP
    conda:
        "../envs/input.yml"
    script:
        "../scripts/handle_provided_files.py"

rule fetch_genome:
    output:
        f"{RESOURCES}/genomes/{{genome}}.fa.gz"
    log:
        f"{LOGS}/genomes/{{genome}}.log"
    benchmark:
        f"{BENCHMARKS}/genomes/{{genome}}.benchmark.txt"
    params:
        path = lambda wildcards: get_genome_path(wildcards.genome, pipeline_config)
    resources:
        tmpdir=TEMP
    conda:
        "../envs/input.yml"
    script:
        "../scripts/fetch_genome.py"
