RESULTS = config['results_path']
RESOURCES = config['resources_path']
LOGS = config['logs_path']
BENCHMARKS = config['benchmarks_path']
TEMP = config['temp_path']

ruleorder: PicardBuildBamIndex > samtools_index

rule unzip:
    input:
        "{file}.gz"
    output:
        "{file}"
    wildcard_constraints:
        file = r"^([\/A-Za-z0-9])*.(fastq|fq|fa|fasta)$"
    resources:
        tmpdir=TEMP
    shell:
        "gzip -dk {input} -f"


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
        cpus_per_task=1,
        mem_mb=4000
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
        BENCHMARKS + "/picard-MarkDuplicates/{sample}.log"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        picard SortSam -I {input.aligned} -O {output.sorted} --SO coordinate
        picard MarkDuplicates -ASO coordinate -I {output.sorted} -O {output.marked} -M {output.metrics} {params.args}
        """
    
'''  
rule samtools_markdup:
    input:
        f"results/{config['aligner']}/{{sample}}.bam"
    output:
        collate = temp("results/samtools-markdup/{sample}_collate.bam"),
        fixmate = temp("results/samtools-markdup/{sample}_fixmate.bam"),
        sort = temp("results/samtools-markdup/{sample}_sort.bam"),
        marked = "results/samtools-markdup/{sample}.bam"
    conda:
        "../envs/utils.yml"
    params:
        args = config['markdup']['args']
    threads:
        12
    shell:
        """
        mkdir -p results/samtools-markdup
        samtools collate -u -o {output.collate} {input} 
        samtools fixmate -@ {threads} -m {output.collate} {output.fixmate}
        samtools sort -@ {threads} -o {output.sort} {output.fixmate}
        samtools markdup -@ {threads} -O BAM {output.sort} {output.marked} 
        """
'''

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
    #log:
        #LOGS + "/samtools-markdup/{sample}.log"
    benchmark:
        BENCHMARKS + "/samtools-markdup/{sample}.txt"
    resources:
        tmpdir=TEMP
    shell:
        """
        #exec > log 2>&1
        mkdir -p {params.path} 
        rm -rf {resources.tmpdir}/{wildcards.sample}_collate &&  rm -rf {resources.tmpdir}/{wildcards.sample}_sort && rm -rf {resources.tmpdir}/{wildcards.sample}_markdup
        samtools collate -T {resources.tmpdir}/{wildcards.sample}_collate -O -u {input} | samtools fixmate -@ {threads} -m -u - - | samtools sort -T {resources.tmpdir}/{wildcards.sample}_sort -@ {threads} -u - | samtools markdup {params.args} -T {resources.tmpdir}/{wildcards.sample}_markdup -@ {threads} - {output}
        """