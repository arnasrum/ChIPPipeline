ruleorder: get_fastq_PE > get_fastq_SE
ruleorder: concatenate_runs_SE > concatenate_runs_PE


rule get_fastq_PE:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        temp(RESOURCES + "/reads/{accession}_1.fastq"),
        temp(RESOURCES + "/reads/{accession}_2.fastq")
    log:
        "logs/pe/{accession}.log"
    params:
        extra="--skip-technical"
    threads: 6  # defaults to 6
    wrapper:
        "v5.10.0/bio/sra-tools/fasterq-dump"

rule get_fastq_SE:
    output:
        temp(RESOURCES + "/reads/{accession}.fastq")
    log:
        "logs/se/{accession}.log"
    params:
        extra="--skip-technical"
    threads: 6
    wrapper:
        "v5.10.0/bio/sra-tools/fasterq-dump"

rule concatenate_runs_SE:
    input:
        lambda wildcards: concatenate_runs_input(file_info["public"][wildcards.sample.split("_")[0]]["runs"], wildcards.sample)
    output:
        [RESOURCES + "/reads/{sample}.fastq.gz"]
    log:
        LOGS + "/concatenate/{sample}.log"
    resources:
        tmpdir=TEMP
    script:
        "../scripts/concatenate_runs.py"

rule concatenate_runs_PE:
    input:
        lambda wildcards: concatenate_runs_input(file_info["public"][wildcards.sample.split("_")[0]]["runs"], wildcards.sample)
    output:
        [RESOURCES + "/reads/{sample}_1.fastq.gz", RESOURCES + "/reads/{sample}_2.fastq.gz"]
    log:
        LOGS + "/concatenate/{sample}.log"
    resources:
        tmpdir=TEMP
    script:
        "../scripts/concatenate_runs.py"

rule handle_provided_samples_SE:
    input:
        lambda wildcards: symlink_input(wildcards.file_name)["read1"]["path"],
    output:
        RESOURCES + "/reads/{file_name}.fastq.gz"
    log:
        LOGS + "/provided_files/{file_name}.log"
    resources:
        tmpdir=TEMP
    script:
        "../scripts/handle_provided_files.py"

rule handle_provided_samples_PE:
    input:
        lambda wildcards: symlink_input(wildcards.file_name)["read1"]["path"],
        lambda wildcards: symlink_input(wildcards.file_name)["read2"]["path"]
    output:
        [f"{RESOURCES}/reads/{{file_name}}_1.fastq.gz", f"{RESOURCES}/reads/{{file_name}}_2.fastq.gz"]
    log:
        LOGS + "/provided_files/{file_name}.log"
    resources:
        tmpdir=TEMP
    script:
        "../scripts/handle_provided_files.py"

rule fetch_genome:
    output:
        RESOURCES + "/genomes/{genome}.fa.gz"
    conda:
        "../envs/download.yml"
    log:
        LOGS + "/genomes/{genome}.log"
    benchmark:
        BENCHMARKS + "/genomes/{genome}.benchmark.txt"
    params:
        genome = lambda wildcards: wildcards.genome
    resources:
        tmpdir=TEMP
    script:
        "../scripts/fetch_genome.py"
