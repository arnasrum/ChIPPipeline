rule trim:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        [f"{RESULTS}/{config['trimmer']}/{{sample}}_1.fastq.gz", f"{RESULTS}/{config['trimmer']}/{{sample}}_2.fastq.gz"]
        if sfs.is_paired_end() else
        [f"{RESULTS}/{config['trimmer']}/{{sample}}.fastq.gz"]
    conda:
        "../envs/trim.yml" if not config["trimmer"] == "cutadapt" else "../envs/cutadapt.yml"
    params:
        args=config[config['trimmer']]["args"],
    threads:
        int(config["trim_threads"])
    log:
        f"{LOGS}/{config['trimmer']}/{{sample}}.log"
    benchmark:
        repeat(f"{BENCHMARKS}/{config['trimmer']}/{{sample}}.txt",config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads
    script:
        "../scripts/trimmers/trim.py"
