
def alignment_input(sample: str) -> list[str]:
    if str(config["paired_end"]).lower() == "true":
        reads = [f"results/{config['trimmer']}/{sample}_1.fastq", f"results/{config['trimmer']}/{sample}_2.fastq"]
    else:
        reads = [f"results/{config['trimmer']}/{sample}.fastq"]
    return reads


rule buildBowtie2Index:
    input:
        f"resources/genomes/{config['genome']}.fa.gz"
    output:
        expand("results/bowtie2-build/{genome}.{extension}", 
            extension=["1.bt2", "2.bt2", "3.bt2", "4.bt2"], 
            genome=[config["genome"]])
    conda:
        "../envs/align.yml"
    params:
        genome = config["genome"],
        args = config["bowtie2-build"]["args"]
    threads:
        int(config["bowtie2-build"]["threads"])
    log:
        f"logs/bowtie2-build/{config['genome']}.log"
    shell:
        '''
        bowtie2-build --threads {threads} {params.args} {input} results/bowtie2-build/{params.genome} 2>&1 {log}
        '''

rule bowtie2:
    input:
        expand("results/bowtie2-build/{genome}.{extension}", 
            extension=["1.bt2", "2.bt2", "3.bt2", "4.bt2"], 
            genome=config["genome"]
        ),
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        f"results/bowtie2/{{sample}}.bam"
    conda:
        "../envs/align.yml"
    params:
        args = config["bowtie2"]["args"],
        extra = config["bowtie2"]["extra"],
        genome = config["genome"],
        paired_end = config["paired_end"]
    threads:
        8
    log:
        "logs/bowtie2/{sample}.log"
    shell:
        '''
        exec > {log} 2>&1
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            inputOptions=''; i=1
            for file in {input.reads}; do inputOptions+="-$i $file "; i=$((i+1)); done
        else
            inputOptions='-q {input.reads}'
        fi
        bowtie2 --mm --threads {threads} -x results/bowtie2-build/{params.genome} $inputOptions {params.args} {params.extra} | samtools view -b -o {output}
        '''

rule buildBWAIndex:
    input:
        f"resources/genomes/{config['genome']}.fa.gz"
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
    log:
        f"logs/bwa-index/{config['genome']}.log"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p results/bwa-index
        bwa index {params.args} {input} -p results/bwa-index/{params.genome}
        """

rule bwa:
    input:
        reads = lambda wildcards: alignment_input(wildcards.sample),
        genomeIndex = expand("results/bwa-index/{genome}.{ext}",
            genome=config["genome"], 
            ext=["amb", "ann", "pac", "sa", "bwt"]
        ) 
    output:
        "results/bwa/{sample}.bam"
    conda:
        "../envs/align.yml"
    params:
        genome = config["genome"],
        args = config["bwa"]["args"],
        extra = config["bwa"]["extra"]
    threads:
        8
    log:
        "logs/bwa-mem/{sample}.log"
    shell:
        """
        exec > {log} 2>&1
        bwa mem -t {threads} {params.args} {params.extra} results/bwa-index/{params.genome} {input.reads} | samtools view -b - > {output}
        """


rule buildStarIndex:
    input:
        f"resources/genomes/{config['genome']}.fa"
    output:
        multiext("results/star-index/SA", "", "index")
    params:
        pathToGenome = f"resources/genomes/{config['genome']}.fa",
        args = config["STAR"]["args"]
    conda:
        "../envs/align.yml"
    threads:
        config["STAR"]["threads"]
    log:
        f"logs/star-index/{config['genome']}.log"
    shell:
        """
        exec > {log} 2>&1
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir results/star-index --genomeFastaFiles {params.pathToGenome} {params.args}
        """ 

rule STAR:
    input:
        expand("results/star-index/SA{index}", index=["", "index"]),
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        "results/STAR/{sample}.bam"
    conda:
        "../envs/align.yml"
    threads:
        config["STAR"]["threads"]
    params:
        args = config["STAR"]["args"],
        extra = config["STAR"]["extra"]
    log:
        "logs/STAR/{sample}.log"
    shell:
        """
        exec > {log} 2>&1
        STAR --readFilesType Fastx --runThreadN {threads} --genomeDir results/star-index --readFilesIn {input.reads} {params.args} {params.extra} --outFileNamePrefix results/STAR/{wildcards.sample} --outStd SAM | samtools view -b -o {output}
        """
