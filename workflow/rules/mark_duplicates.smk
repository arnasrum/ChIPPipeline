rule picard_MarkDuplicates:
    input:
        aligned = RESULTS + "/" + aligner_name + "/{sample}.bam",
        aligned_index = RESULTS + "/" + aligner_name + "/{sample}.bam.bai",
    output:
        sorted = temp(RESULTS + "/MarkDuplicates/{sample}_sorted.bam"),
        marked = RESULTS + "/MarkDuplicates/{sample}.bam",
        metrics = RESULTS + "/MarkDuplicates/{sample}.metrics.txt"
    params:
        args = config['MarkDuplicates']['args']
    conda:
        "../envs/utils.yml"
    log:
        LOGS + "/MarkDuplicates/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/MarkDuplicates/{sample}.log", config["benchmark_repeat_duplicate"])
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        picard SortSam --TMP_DIR {resources.tmpdir} -I {input.aligned} -O {output.sorted} --SO coordinate
        picard MarkDuplicates -ASO coordinate --TMP_DIR {resources.tmpdir} -I {output.sorted} -O {output.marked} -M {output.metrics} {params.args}
        """

rule samtools_markdup:
    input:
        RESULTS + "/" + aligner_name + "/{sample}.bam"
    output:
        RESULTS + "/markdup/{sample}.bam",
    conda:
        "../envs/utils.yml"
    params:
        args = config['markdup']['args'],
        path = f"{RESULTS}/markdup"
    threads:
        int(config["markdup"]["threads"])
    log:
        LOGS + "/markdup/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/markdup/{sample}.txt", config["benchmark_repeat_duplicate"])
    resources:
        tmpdir=TEMP,
        cpus_per_thread=lambda wildcards, threads: threads
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.path} 
        uuid=$(python3 -c "import uuid; print(uuid.uuid4())")
        samtools collate -O -T {resources.tmpdir}/${{uuid}}_collate {input} \
            | samtools fixmate -@ {threads} -m - - \
            | samtools sort -T {resources.tmpdir}/${{uuid}}_sort -@ {threads} - -u \
            | samtools markdup {params.args} -T {resources.tmpdir}/${{uuid}}_markdup -@ {threads} - {output}
        """