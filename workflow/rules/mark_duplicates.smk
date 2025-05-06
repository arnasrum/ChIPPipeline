rule picard_MarkDuplicates:
    input:
        aligned = f"{RESULTS}/{config['aligner']}/{{sample}}.bam",
    output:
        sorted = temp(f"{RESULTS}/MarkDuplicates/{{sample}}_sorted.bam"),
        marked = f"{RESULTS}/MarkDuplicates/{{sample}}.bam",
        index = f"{RESULTS}/MarkDuplicates/{{sample}}.bam.bai",
        metrics = f"{RESULTS}/MarkDuplicates/{{sample}}.metrics.txt"
    params:
        args = config['MarkDuplicates']['args']
    conda:
        "../envs/utils.yml"
    log:
        f"{LOGS}/MarkDuplicates/{{sample}}.log"
    benchmark:
        repeat(f"{BENCHMARKS}/MarkDuplicates/{{sample}}.txt", config["benchmark_repeat_duplicate"])
    resources:
        tmpdir=TEMP,
        cpus_per_thread= lambda wildcards,threads: threads,
        mem_mb= lambda wildcards,attempt: config["MarkDuplicates"]["mem_mb"] * attempt,
        runtime= lambda wildcards,attempt: config["MarkDuplicates"]["runtime"] * attempt,
    shell:
        """
        exec > {log} 2>&1
        picard SortSam --TMP_DIR {resources.tmpdir} -I {input.aligned} -O {output.sorted} --SO coordinate
        picard MarkDuplicates -ASO coordinate --TMP_DIR {resources.tmpdir} -I {output.sorted} -O {output.marked} -M {output.metrics} {params.args}
        picard BuildBamIndex -I {output.marked} -O {output.index}
        """

rule samtools_markdup:
    input:
        f"{RESULTS}/{config['aligner']}/{{sample}}.bam"
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
        repeat(f"{BENCHMARKS}/markdup/{{sample}}.txt", config["benchmark_repeat_duplicate"])
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