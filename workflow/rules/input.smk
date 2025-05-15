ruleorder: fasterq_dump_PE > fasterq_dump_SE
ruleorder: concatenate_runs_PE > concatenate_runs_SE

rule fasterq_dump_PE:
    output:
        temp(f"{RESOURCES}/reads/{{accession}}_1.fastq"),
        temp(f"{RESOURCES}/reads/{{accession}}_2.fastq")
    log:
        f"{LOGS}/fasterq-dump_pe/{{accession}}.log"
    params:
        args=str(config["fasterq-dump"]["args"]),
        outdir=f"{RESOURCES}/reads"
    threads:
        int(config["fasterq-dump"]["threads"])
    resources:
        tmpdir = TEMP,
        cpus_per_task = int(config["fasterq-dump"]["threads"]),
        runtime = int(config["fasterq-dump"]["runtime"]),
        mem_mb = int(config["fasterq-dump"]["mem_mb"])
    conda:
        "../envs/input.yml"
    shell:
        """
            fasterq-dump {wildcards.accession} {params.args} -e {threads} -t {resources.tmpdir} -O {params.outdir} -F fastq
        """

rule fasterq_dump_SE:
    output:
        temp(f"{RESOURCES}/reads/{{accession}}.fastq")
    log:
        f"{LOGS}/fasterq-dump_se/{{accession}}.log"
    params:
        args=lambda wildcards, resources: str(config["fasterq-dump"]["args"]),
        outdir=f"{RESOURCES}/reads"
    threads:
        int(config["fasterq-dump"]["threads"])
    resources:
        tmpdir = TEMP,
        cpus_per_task = int(config["fasterq-dump"]["threads"]),
        runtime = int(config["fasterq-dump"]["runtime"]),
        mem_mb = int(config["fasterq-dump"]["mem_mb"])
    shell:
        """
        fasterq-dump {wildcards.accession} {params.args} -e {threads} -t {resources.tmpdir} -O {params.outdir} -F fastq
        """

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
        "../scripts/tools/concatenate_runs.py"

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
        "../scripts/tools/concatenate_runs.py"

rule handle_provided_samples_SE:
    input:
        lambda wildcards: handle_provided_input(wildcards.file_name, pipeline_config),
    output:
        f"{RESOURCES}/reads/{{file_name}}.fastq.gz"
    log:
        f"{LOGS}/provided_files/{{file_name}}.log"
    resources:
        tmpdir=TEMP
    conda:
        "../envs/input.yml"
    params:
        outdir = f"{RESOURCES}/reads"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        if [[ {input} = *.gz]]; then
            ln -srf {input} {output}
        else
            gzip -kfc {input} > {output} 
        fi
        """

rule handle_provided_samples_PE:
    input:
        lambda wildcards: handle_provided_input(wildcards.file_name, pipeline_config),
    output:
        f"{RESOURCES}/reads/{{file_name}}_1.fastq.gz",
        f"{RESOURCES}/reads/{{file_name}}_2.fastq.gz"
    log:
        f"{LOGS}/provided_files/{{file_name}}.log"
    params:
        outdir = f"{RESOURCES}/reads"
    resources:
        tmpdir=TEMP
    conda:
        "../envs/input.yml"
    shell:
        """
            exec > {log} 2>&1
            mkdir -p {params.outdir}
            input="{input}"  
            output="{output}"  
            INPUT_ARRAY=( $input )
            OUTPUT_ARRAY=( $output )
            OUTPUT_DIR="{params.outdir}"
            for index in "${{!INPUT_ARRAY[@]}}"; do
                input_file="${{INPUT_ARRAY[index]}}"
                output_file="${{OUTPUT_ARRAY[index]}}"
                if [[ "$input_file" = *.gz ]]; then
                    ln -srf "$input_file" "$output_file"
                else
                    gzip -kfc "$input_file" > "$output_file"
                fi
            done
        """

rule fetch_genome:
    output:
        f"{RESOURCES}/genomes/{{genome}}.fa.gz"
    log:
        f"{LOGS}/genomes/{{genome}}.log"
    benchmark:
        f"{BENCHMARKS}/genomes/{{genome}}.benchmark.txt"
    params:
        path = lambda wildcards: get_genome_path(wildcards.genome, pipeline_config),
    resources:
        tmpdir=TEMP
    conda:
        "../envs/input.yml"
    shell:
        """
            exec > {log} 2>&1
            if [[ "{params.path}" = "None" ]]; then
                echo "Defined genome is not a file; trying to download"
                rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz {output}
            else
                if [[ "{params.path}" = *.gz ]]; then
                    echo "gzipped genome detected; symlinking"
                    ln -sr {params.path} {output}
                else
                    echo "non-gzipped genome detected; gzipping"
                    gzip -kfc {params.path} > {output}
                fi
            fi
        """

