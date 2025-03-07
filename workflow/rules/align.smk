rule build_index:
    input:
        build_index_input()
    output:
        index_files
    conda:
        "../envs/align.yml"
    params:
        args=config["index_args"],
    threads:
        int(config["indexing_threads"])
    log:
        f"{LOGS}/{aligner.get_name()}_index/{genome}.log"
    benchmark:
        f"{BENCHMARKS}/{aligner.get_name()}_index/{genome}.txt"
    resources:
        tmpdir=TEMP,
    script:
        "../scripts/aligners/build_index.py"

rule align:
    input:
        index = index_files,
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        f"{RESULTS}/{aligner_name}/{{sample}}.bam"
    conda:
        "../envs/align.yml"
    params:
        args=config[aligner_name]["args"],
    threads:
        int(config["aligning_threads"])
    log:
        f"{LOGS}/{aligner_name}/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/{aligner_name}/{{sample}}.txt"
    resources:
        tmpdir=TEMP,
    script:
        "../scripts/aligners/align.py"
