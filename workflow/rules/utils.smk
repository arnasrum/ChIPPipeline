RESULTS = config['results_path']
RESOURCES = config['resources_path']
LOGS = config['logs_path']
BENCHMARKS = config['benchmarks_path']
TEMP = config['temp_path']

rule unzip:
    input:
        "{file}.gz"
    output:
        "{file}"
    wildcard_constraints:
        file = r"^([\/A-Za-z0-9])*.(fastq|fq|fa|fasta).gz"
    resources:
        tmpdir=TEMP
    shell:
        "gzip -dk {input} -f"


rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    threads:
        int(config["samtools-index"]["threads"])
    conda:
        "../envs/utils.yml"
    resources:
        tmpdir=TEMP
    shell:
        """
        samtools index -@ {threads} {input} -o {output}
        """

rule picardCreateGenomeSequenceDictionary:
    input:
        RESOURCES + "/genomes/{genome}.fa"
    output:
        RESOURCES + "/genomes/{genome}.dict"
    conda:
        "../envs/utils.yml"
    resources:
        tmpdir=TEMP
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output}"

rule picard_MarkDuplicates:
    input:
        aligned = f"{RESULTS}/{config['aligner']}/{{sample}}.bam",
        aligned_index= f"{RESULTS}/{config['aligner']}/{{sample}}.bam.bai",
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
        BENCHMARKS + "/picard-MarkDuplicates/{sample}.log"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        picard SortSam -I {input.aligned} -O {output.sorted} --SO coordinate
        picard MarkDuplicates -I {output.sorted} -O {output.marked} -M {output.metrics} {params.args}
        """

rule samtools_markdup:
    input:
        f"{RESULTS}/{config['aligner']}/{{sample}}.bam"
    output:
        RESULTS + "/samtools-markdup/{sample}.bam"
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
        BENCHMARKS + "/samtools-markdup/{sample}.txt"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.path} 
        samtools collate -T {resources.tmpdir}/{wildcards.sample}_collate -@ {threads} -O -u {input} | samtools fixmate -@ {threads} -m -u - - | samtools sort -T {resources.tmpdir}/{wildcards.sample}_sort -@ {threads} -u - | samtools markdup {params.args} -T {resources.tmpdir}/{wildcards.sample}_markdup -@ {threads} - {output}
        """