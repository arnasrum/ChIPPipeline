ruleorder: fastq_dump_SE > fastq_dump_PE
ruleorder: concatenate_runs_SE > concatenate_runs_PE

rule fastq_dump_SE:
    output:
        temp(RESOURCES + "/reads/{srr}.fastq")
    params:
        path = f"{RESOURCES}/reads"
    conda:
        "../envs/download.yml"
    wildcard_constraints:
        srr = r"SRR[0-9]*"
    log:
        LOGS + "/fastq-dump/{srr}_SE.log"
    benchmark:
        BENCHMARKS + "/fastq-dump/{srr}_SE.log"
    resources:
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        fastq-dump -O {params.path} {wildcards.srr}
        '''

rule fastq_dump_PE:
    output:
        temp(RESOURCES + "/reads/{srr}_1.fastq"),
        temp(RESOURCES + "/reads/{srr}_2.fastq")
    wildcard_constraints:
        srr = r"SRR[0-9]*"
    params:
        path = f"{RESOURCES}/reads"
    conda:
        "../envs/download.yml"
    log:
        LOGS + "/fastq-dump/{srr}_PE.log"
    benchmark:
        BENCHMARKS + "/fastq-dump/{srr}_PE.log"
    resources:
        tmpdir=TEMP
    threads: 6
    shell:
        '''
        exec > {log} 2>&1
        fasterq-dump -t {resources.tmpdir} -e {threads} -O {params.path} --split-files {wildcards.srr}
        '''

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
