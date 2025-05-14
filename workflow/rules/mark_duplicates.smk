rule picard_mark_duplicates:
    input:
        aligned = f"{RESULTS}/{pipeline_config.get_config_option('aligner')}/{{sample}}.bam",
    output:
        sorted = temp(f"{RESULTS}/mark_duplicates/{{sample}}_sorted.bam"),
        marked = f"{RESULTS}/mark_duplicates/{{sample}}.bam",
        index = f"{RESULTS}/mark_duplicates/{{sample}}.bam.bai",
        metrics = f"{RESULTS}/mark_duplicates/{{sample}}.metrics.txt"
    params:
        args = config['mark_duplicates']['args']
    conda:
        "../envs/utils.yml"
    log:
        f"{LOGS}/mark_duplicates/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/mark_duplicates/{{sample}}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_thread= lambda wildcards,threads: threads,
        mem_mb= lambda wildcards,attempt: config["mark_duplicates"]["mem_mb"] * attempt,
        runtime= lambda wildcards,attempt: config["mark_duplicates"]["runtime"] * attempt,
    shell:
        """
        exec > {log} 2>&1
        picard SortSam --TMP_DIR {resources.tmpdir} -I {input.aligned} -O {output.sorted} --SO coordinate
        picard MarkDuplicates -ASO coordinate --TMP_DIR {resources.tmpdir} -I {output.sorted} -O {output.marked} -M {output.metrics} {params.args}
        picard BuildBamIndex -I {output.marked} -O {output.index}
        """

rule samtools_markdup:
    input:
        f"{RESULTS}/{pipeline_config.get_config_option('aligner')}/{{sample}}.bam"
    output:
        marked = f"{RESULTS}/markdup/{{sample}}.bam",
        index = f"{RESULTS}/markdup/{{sample}}.bam.bai",
        stats = f"{RESULTS}/markdup/{{sample}}_stats.txt"
    conda:
        "../envs/utils.yml"
    params:
        args = config['markdup']['args'],
        path = f"{RESULTS}/markdup"
    threads:
        int(config["markdup"]["threads"])
    log:
        f"{LOGS}/markdup/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/markdup/{{sample}}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_thread=lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: config["markdup"]["mem_mb"] * attempt,
        runtime = lambda wildcards,attempt: config["markdup"]["runtime"] * attempt,
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.path} 
        uuid=$(python3 -c "import uuid; print(uuid.uuid4())")
        samtools collate -@ {threads} -u -O -T {resources.tmpdir}/${{uuid}}_collate {input} \
            | samtools fixmate -@ {threads} -u -m - - \
            | samtools sort -T {resources.tmpdir}/${{uuid}}_sort -@ {threads} - -u \
            | samtools markdup {params.args} -f {output.stats} -T {resources.tmpdir}/${{uuid}}_markdup -@ {threads} - {output.marked}
        samtools index -@ {threads} {output.marked} {output.index}
        """