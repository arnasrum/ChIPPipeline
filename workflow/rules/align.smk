trimmer = config["trimmer"]
aligner = config["aligner"]


def alignmentInput(sample: str) -> list[str]: 
    reads = [f"results/{config["trimmer"]}/{sample}_1.fastq"]
    if config["paired_end"]:
        reads.append(f"results/{config["trimmer"]}/{sample}_2.fastq") 
    return reads


rule buildBowtie2Index:
    input:
        f"resources/genomes/{config["genome"]}.fa.gz"
    output:
        expand("results/bowtie2-build/{genome}.{extension}", 
            extension=["1.bt2", "2.bt2", "3.bt2", "4.bt2"], 
            genome=[config["genome"]])
    conda:
        "../envs/align.yml"
    params:
        genome = config["genome"],
        args = config["bowtie2-build"]["args"]
    shell:
        '''
        bowtie2-build {params.args} {input} results/bowtie2-build/{params.genome}
        '''

rule bowtie2:
    input:
        expand("results/bowtie2-build/{genome}.{extension}", 
            extension=["1.bt2", "2.bt2", "3.bt2", "4.bt2"], 
            genome=config["genome"]
        ),
        reads = lambda wildcards: alignmentInput(wildcards.sample)
    output:
        outputFile = f"results/bowtie2/{{sample}}.sam"
    conda:
        "../envs/align.yml"
    params:
        args = config["bowtie2"]["args"],
        extra = config["bowtie2"]["extra"],
        genome = config["genome"],
        paired_end = config["paired_end"],
    log:
        "logs/bowtie2/{sample}.log"
    shell:
        '''
        if [[ params.paired_end ]]; then
            bowtie2 -x results/bowtie2-build/{params.genome} -1 {input.reads[0]} -2 {input.reads[1]} -S {output} {params.args} {params.extra}
        else
            bowtie2 -x results/bowtie2-build/{params.genome} -1 {input.reads} -S {output} {params.args} {params.extra}
        fi
        '''

rule buildBWAIndex:
    input:
        f"resources/genomes/{config["genome"]}.fa.gz"
    output: 
        expand("results/bwa-index/{genome}.{ext}", 
            genome=config["genome"], 
            ext=["amb", "ann", "pac", "sa", "bwt"]
        ) 
    conda:
        "../envs/align.yml"
    params:
        genome = config["genome"],
        args = config["bwa-index"]["args"]
    shell:
        """
        mkdir -p results/bwa-index
        bwa index {params.args} {input} -p results/bwa-index/{params.genome}
        """

rule bwa:
    input:
        reads = lambda wildcards: alignmentInput(wildcards.sample),
        genomeIndex = expand("results/bwa-index/{genome}.{ext}",
            genome=config["genome"], 
            ext=["amb", "ann", "pac", "sa", "bwt"]
        ) 
    output:
        "results/bwa/{sample}.sam"
    conda:
        "../envs/align.yml"
    params:
        args = config["bwa"]["args"],
        extra = config["bwa"]["extra"]
    shell:
        """
        bwa mem {params.args} results/bwa-index/mm39 {input.reads} > {output} {params.extra}
        """


rule buildStarIndex:
    input:
        f"resources/genomes/{config["genome"]}.fa"
    output:
        expand("results/star-index/SA{index}", index=["", "index"])
    params:
        pathToGenome = f"resources/genomes/{config["genome"]}.fa"
    conda:
        "../envs/align.yml"
    threads:
        config["STAR"]["threads"]
    params:
        args = config["STAR"]["args"]
    shell:
        """
        mkdir -p results/starIndex
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir results/star-index --genomeFastaFiles {params.pathToGenome} {params.args}
        """ 

rule STAR:
    input:
        expand("results/star-index/SA{index}", index=["", "index"]),
        reads = lambda wildcards: alignmentInput(wildcards.sample)
    output:
        "results/STAR/{sample}.sam"
    conda:
        "../envs/align.yml"
    threads:
        config["STAR"]["threads"]
    params:
        args = config["STAR"]["args"],
        extra = config["STAR"]["extra"]
    shell:
        """
        STAR --runThreadN {threads} --genomeDir results/star-index --outFileNamePrefix results/STAR/{wildcards.sample} --readFilesIn {input.reads} {params.args} {params.extra}
        mv results/STAR/{wildcards.sample}Aligned.out.sam results/STAR/{wildcards.sample}.sam
        """

rule filterReads:
    input:
        f"results/{config["aligner"]}/{{id}}.sam"
    output:
        f"results/samtools/{{id}}_filtered.bam"
    params:
        args = config["samtools"],
        aligner = config["aligner"]
    log:
        "logs/samtools/{id}.log"
    shell:
        '''
        if [ ! -d "results/samtools" ]
        then
            mkdir "results/samtools"
        fi
        samtools view {params.args} -b -o {output} {input} 
        '''