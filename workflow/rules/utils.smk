
ruleorder: PicardBuildBamIndex > samtools_index

rule unzip:
    input:
        "{file}.gz"
    output:
        "{file}"
    wildcard_constraints:
        file = r"^(.*).(fastq|fq|fa|fasta)$"
    resources:
        tmpdir=TEMP
    params:
        args = config["gzip"]["args"]
    shell:
        "gzip {params.args} -dkf {input}"

rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bai"
    threads:
        int(config["samtools-index"]["threads"])
    log:
        LOGS + "/samtools-index/{sample}.log"
    conda:
        "../envs/utils.yml"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        samtools index -@ {threads} - -o {output}
        """


rule PicardBuildBamIndex:
    input:
        "{sample}.bam"
    output:
        "{sample}.bai"
    log:
        LOGS + "/BuildBamIndex/{sample}.log"
    conda:
        "../envs/utils.yml"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        picard BuildBamIndex -I {input} -O {output}
        """


rule picardCreateGenomeSequenceDictionary:
    input:
        RESOURCES + "/genomes/{genome}.fa"
    output:
        RESOURCES + "/genomes/{genome}.dict"
    conda:
        "../envs/utils.yml"
    resources:
        tmpdir=TEMP,
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output}"

rule picard_MarkDuplicates:
    input:
        aligned = RESULTS + "/" + config['aligner'] + "/{sample}.bam",
        aligned_index = RESULTS + "/" + config['aligner'] + "/{sample}.bam.bai",
    output:
        sorted = temp(RESULTS + "/picard-MarkDuplicates/{sample}_sorted.bam"),
        marked = RESULTS + "/picard-MarkDuplicates/{sample}.bam",
        metrics = RESULTS + "/picard-MarkDuplicates/{sample}.metrics.txt"
    params:
        paired_end = config['paired_end'],
        genome = config['genome'],
        aligner = config['aligner'],
        args = config['MarkDuplicates']['args']
    conda:
        "../envs/utils.yml"
    log:
        LOGS + "/picard-MarkDuplicates/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/picard-MarkDuplicates/{sample}.log", config["benchmark_repeat_duplicate"])
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        picard SortSam -I {input.aligned} -O {output.sorted} --SO coordinate
        picard MarkDuplicates -ASO coordinate -I {output.sorted} -O {output.marked} -M {output.metrics} {params.args}
        """

rule samtools_markdup:
    input:
        RESULTS + "/" + config['aligner'] + "/{sample}.bam"
    output:
        RESULTS + "/samtools-markdup/{sample}.bam",
    conda:
        "../envs/utils.yml"
    params:
        args = config['markdup']['args'],
        path = f"{RESULTS}/samtools-markdup"
    threads:
        int(config["markdup"]["threads"])
    log:
        LOGS + "/samtools-markdup/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/samtools-markdup/{sample}.txt", config["benchmark_repeat_duplicate"])
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.path} 
        uuid=$(python3 -c "import uuid; print(uuid.uuid4())")
        samtools collate -T {resources.tmpdir}/{{uuid}}_collate -O -u {input} \
            | samtools fixmate -@ {threads} -m -u - - \
            | samtools sort -T {resources.tmpdir}/{{uuid}}_sort -@ {threads} -u - \
            | samtools markdup {params.args} -T {resources.tmpdir}/{{uuid}}_markdup -@ {threads} - {output}
        """