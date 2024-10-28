trimmer = config["trimmer"]
aligner = config["aligner"]


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

rule bowtie2_pe:
    input:
        expand("results/bowtie2-build/{genome}.{extension}", 
            extension=["1.bt2", "2.bt2", "3.bt2", "4.bt2"], 
            genome=config["genome"]
        ),
        input1 = f"results/{trimmer}/{{id}}_1.fastq", 
        input2 = f"results/{trimmer}/{{id}}_2.fastq"
    output:
        outputFile = f"results/bowtie2/{{id}}.sam"
    conda:
        "../envs/align.yml"
    params:
        args = config["bowtie2"]["args"],
        genome = config["genome"],
        outputArgs = f"-S results/bowtie2/{{id}}.sam"
    log:
        "logs/bowtie2/{id}.log"
    shell:
        '''
        bowtie2 -x results/bowtie2-build/{params.genome} -1 {input.input1} -2 {input.input2} -S {output}
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

rule bwa_pe:
    input:
        read1 = f"results/{config["trimmer"]}/{{id}}_1.fastq",
        read2 = f"results/{config["trimmer"]}/{{id}}_2.fastq",
        genomeIndex = expand("results/bwa-index/{genome}.{ext}", 
            genome=config["genome"], 
            ext=["amb", "ann", "pac", "sa", "bwt"]
        ) 
    output:
        "results/bwa/{id}.sam"
    conda:
        "../envs/align.yml"
    params:
        args = config["bwa"]["args"]
    shell:
        """
        bwa mem {params.args} results/bwa-index/mm39 {input.read1} {input.read2} > {output}
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

rule STAR_pe:
    input:
        expand("results/star-index/SA{index}", index=["", "index"]),
        read1 = f"results/{config["trimmer"]}/{{sample}}_1.fastq",
        read2 = f"results/{config["trimmer"]}/{{sample}}_2.fastq"
    output:
        "results/STAR/{sample}.sam"
    conda:
        "../envs/align.yml"
    threads:
        config["STAR"]["threads"]
    params:
        args = config["STAR"]["args"]
    shell:
        """
        STAR --runThreadN {threads} --genomeDir results/star-index --outFileNamePrefix results/STAR/{wildcards.sample} --readFilesIn {input.read1} {input.read2} {params.args}
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