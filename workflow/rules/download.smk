ruleorder: symlink_reference_genome > download_reference_genome
localrules: download_reference_genome, symlink_reference_genome

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



read_extensions = ["_1", "_2"] if sfs.is_paired_end() else [""]
rule concatenate_runs:
    input:
        lambda wildcards: concatenate_runs_input(file_info["public"][wildcards.sample.split("_")[0]]["runs"], read_extensions)
    output:
        [RESOURCES + "/reads/{sample}_1.fastq.gz", RESOURCES + "/reads/{sample}_2.fastq.gz"]
        if sfs.is_paired_end() else
        [RESOURCES + "/reads/{sample}.fastq.gz"]
    log:
        LOGS + "/concatenate/{sample}.log"
    resources:
        tmpdir=TEMP
    script:
        "../scripts/concatenate_runs.py"

rule handle_provided_files:
    input:
        lambda wildcards: symlink_input(wildcards.file_name)["read1"]["path"],
        lambda wildcards: symlink_input(wildcards.file_name)["read2"]["path"] if sfs.is_paired_end() else ""
    output:
        [f"{RESOURCES}/reads/{{file_name}}_1.fastq.gz", f"{RESOURCES}/reads/{{file_name}}_2.fastq.gz"]
        if sfs.is_paired_end() else
        RESOURCES + "/reads/{file_name}.fastq.gz"
    log:
        LOGS + "/provided_files/{file_name}.log"
    resources:
        tmpdir=TEMP
    script:
        "../scripts/handle_provided_files.py"

rule symlink_reference_genome:
    input:
        reference_genome_input()
    output:
        RESOURCES + "/genomes/{genome}.fa.gz"
    params:
        gzip_args = config["gzip"]["args"]
    log:
        LOGS + "/ln/{genome}.log"
    benchmark:
        BENCHMARKS + "/ln/{genome}.benchmark.txt"
    resources:
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        if [[ {input} == *.gz ]]; then
            echo "The string ends with .gz, skipping compression."
            ln -sr {input} {output} 
        else
            echo "The string does not end with .gz, compressing genome."
            gzip {params.gzip_args} -kfc {input} > {output}
        fi 
        '''

rule download_reference_genome:
    output:
        RESOURCES + "/genomes/{genome}.fa.gz"
    conda:
        "../envs/download.yml"
    log:
        LOGS + "/rsync/{genome}.log"
    benchmark:
        BENCHMARKS + "/rsync/{genome}.benchmark.txt"
    resources:
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz {output}
        '''
