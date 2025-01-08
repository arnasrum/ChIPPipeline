import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import SampleFileScripts

sfs = SampleFileScripts(config)

RESULTS: str = config['results_path']
RESOURCES: str = config['resources_path']
LOGS: str = config['logs_path']
BENCHMARKS: str = config['benchmarks_path']
TEMP: str = config['temp_path']


def alignment_input(sample: str) -> list[str]:
    file_info = SampleFileScripts.get_file_info(config["json_path"])
    if sfs.is_paired_end():
        if sample in [*map(lambda accession: file_info['public'][accession]['file_name'], file_info["public"].keys())]:
            reads = [f"{RESULTS}/{config['trimmer']}/{sample}_1.fastq", f"{RESULTS}/{config['trimmer']}/{sample}_2.fastq"]
        else:
            reads = [f"{RESULTS}/{config['trimmer']}/{file_info['provided'][sample]['file_name']}_1.fastq",
                     f"{RESULTS}/{config['trimmer']}/{file_info['provided'][sample]['file_name']}_2.fastq"]
    else:
        if sample in [*map(lambda accession: file_info['public'][accession]['file_name'],file_info["public"].keys())]:
            reads = [f"{RESULTS}/{config['trimmer']}/{sample}.fastq"]
        else:
            reads = [f"{RESULTS}/{config['trimmer']}/{file_info['provided'][sample]['file_name']}.fastq"]
    return reads

rule buildBowtie2Index:
    input:
        f"{RESOURCES}/genomes/{config['genome']}.fa.gz"
    output:
        multiext(f"{RESULTS}/bowtie2-build/{config['genome']}.", "1.bt2", "2.bt2", "3.bt2", "4.bt2"),
    conda:
        "../envs/align.yml"
    params:
        genome = config["genome"],
        args = config["bowtie2-build"]["args"],
        path = RESULTS + "/bowtie2-build"
    threads:
        int(config["bowtie2-build"]["threads"])
    log:
        f"{LOGS}/bowtie2-build/{config['genome']}.log"
    benchmark:
        f"{BENCHMARKS}/bowtie2-build/{config['genome']}.txt"
    resources:
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        bowtie2-build --threads {threads} {params.args} {input} {params.path}/{params.genome}
        '''

rule bowtie2:
    input:
        multiext(f"{RESULTS}/bowtie2-build/{config['genome']}.", "1.bt2", "2.bt2", "3.bt2", "4.bt2"),
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        temp(RESULTS + "/bowtie2/{sample}.bam")
    conda:
        "../envs/align.yml"
    params:
        args = config["bowtie2"]["args"],
        extra = config["bowtie2"]["extra"],
        genome = config["genome"],
        paired_end = config["paired_end"],
        index_path = f"{RESULTS}/bowtie2-build/{config['genome']}"
    threads:
        int(config["bowtie2"]["threads"])
    log:
        LOGS + "/bowtie2/{sample}.log"
    benchmark:
        BENCHMARKS + "/bowtie2/{sample}.txt"
    resources:
        tmpdir=TEMP
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
        bowtie2 --mm --threads {threads} -x {params.index_path} $inputOptions {params.args} {params.extra} | samtools view -b -o {output}
        '''

rule buildBWAIndex:
    input:
        f"{RESOURCES}/genomes/{config['genome']}.fa.gz"
    output: 
        multiext(f"{RESULTS}/bwa-index/{config['genome']}.", "amb", "ann", "pac", "sa", "bwt")
    conda:
        "../envs/align.yml"
    params:
        index_path = f"{RESULTS}/bwa-index/{config['genome']}",
        dir = f"{RESULTS}/bwa-index/",
        args = config["bwa-index"]["args"]
    log:
        f"{LOGS}/bwa-index/{config['genome']}.log"
    benchmark:
        f"{BENCHMARKS}/bwa-index/{config['genome']}.txt"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.dir} 
        bwa index {params.args} {input} -p {params.index_path} 
        """

rule bwa:
    input:
        reads = lambda wildcards: alignment_input(wildcards.sample),
        genomeIndex = multiext(f"{RESULTS}/bwa-index/{config['genome']}.", "amb", "ann", "pac", "sa", "bwt")
    output:
        temp(RESULTS + "/bwa/{sample}.bam")
    conda:
        "../envs/align.yml"
    params:
        genome = config["genome"],
        args = config["bwa"]["args"],
        extra = config["bwa"]["extra"],
        index_path= f"{RESULTS}/bwa-index/config['genome']",
    threads:
        int(config['bwa']['threads'])
    log:
        LOGS + "/bwa-mem/{sample}.log"
    benchmark:
        BENCHMARKS + "/bwa-mem/{sample}.txt"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        bwa mem -t {threads} {params.args} {params.extra} {params.index_path} {input.reads} | samtools view -b - > {output}
        """


rule buildStarIndex:
    input:
        f"{RESOURCES}/genomes/{config['genome']}.fa"
    output:
        multiext(RESULTS + "/star-index/SA", "", "index")
    params:
        genome_path = f"{RESOURCES}/genomes/{config['genome']}.fa",
        result_path = f"{RESULTS}/star-index",
        args = config["STAR"]["args"]
    conda:
        "../envs/align.yml"
    threads:
        config["STAR"]["threads"]
    log:
        f"{LOGS}/star-index/{config['genome']}.log"
    benchmark:
        f"{BENCHMARKS}/star-index/{config['genome']}.txt"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.result_path} --genomeFastaFiles {params.genome_path} {params.args}
        """ 

rule STAR:
    input:
        multiext(RESULTS + "/star-index/SA", "", "index"),
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        temp(RESULTS + "/STAR/{sample}.bam")
    conda:
        "../envs/align.yml"
    threads:
        config["STAR"]["threads"]
    params:
        args = config["STAR"]["args"],
        extra = config["STAR"]["extra"],
        index_path = f"{RESULTS}/star-index"
    log:
        LOGS + "/STAR/{sample}.log"
    benchmark:
        BENCHMARKS + "/STAR/{sample}.txt"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        #STAR --readFilesType Fastx --runThreadN {threads} --genomeDir {params.index_path} --readFilesIn {input.reads} {params.args} {params.extra} --outFileNamePrefix results/STAR/{wildcards.sample} --outStd SAM | samtools view -b -o {output}
        STAR --readFilesType Fastx --runThreadN {threads} --genomeDir {params.index_path} --readFilesIn {input.reads} {params.args} {params.extra} --outStd SAM | samtools view -b -o {output}
        """
