genome = sfs.get_genome_code()

rule buildBowtie2Index:
    input:
        f"{RESOURCES}/genomes/{genome}.fa.gz"
    output:
        multiext(f"{RESULTS}/bowtie2-build/{genome}.", "1.bt2", "2.bt2", "3.bt2", "4.bt2"),
    conda:
        "../envs/align.yml"
    params:
        genome = genome,
        args = config["bowtie2-build"]["args"],
        path = RESULTS + "/bowtie2-build"
    threads:
        int(config["bowtie2-build"]["threads"])
    log:
        f"{LOGS}/bowtie2-build/{genome}.log"
    benchmark:
        f"{BENCHMARKS}/bowtie2-build/{genome}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_task=int(config["bowtie2-build"]["threads"]),
        mem_mb=20000
    shell:
        '''
        exec > {log} 2>&1
        bowtie2-build --threads {threads} {params.args} {input} {params.path}/{params.genome}
        '''

rule bowtie2:
    input:
        multiext(f"{RESULTS}/bowtie2-build/{genome}.", "1.bt2", "2.bt2", "3.bt2", "4.bt2"),
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        RESULTS + "/bowtie2/{sample}.bam"
    conda:
        "../envs/align.yml"
    params:
        args = config["bowtie2"]["args"],
        extra = config["bowtie2"]["extra"],
        genome = config["genome"],
        paired_end = config["paired_end"],
        index_path = f"{RESULTS}/bowtie2-build/{genome}"
    threads:
        int(config["bowtie2"]["threads"])
    log:
        LOGS + "/bowtie2/{sample}.log"
    benchmark:
        BENCHMARKS + "/bowtie2/{sample}.txt"
    resources:
        tmpdir=TEMP,
        mem_per_cpu=2096,
        cpus_per_task=int(config["bowtie2"]["threads"])
    shell:
        '''
        exec > {log} 2>&1
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            inputOptions=''; i=1
            for file in {input.reads}; do inputOptions+="-$i $file "; i=$((i+1)); done
        else
            inputOptions='-U {input.reads}'
        fi
        bowtie2 -q --threads {threads} -x {params.index_path} $inputOptions {params.args} {params.extra} | samtools sort -@ {threads} -o {output} -
        '''

rule buildBWAIndex:
    input:
        f"{RESOURCES}/genomes/{genome}.fa.gz"
    output: 
        multiext(f"{RESULTS}/bwa-index/{genome}.", "amb", "ann", "pac", "sa", "bwt")
    conda:
        "../envs/align.yml"
    params:
        index_path = f"{RESULTS}/bwa-index/{genome}",
        dir = f"{RESULTS}/bwa-index/",
        args = config["bwa-index"]["args"]
    log:
        f"{LOGS}/bwa-index/{genome}.log"
    benchmark:
        f"{BENCHMARKS}/bwa-index/{genome}.txt"
    resources:
        tmpdir=TEMP,
        mem_mb=16000
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.dir} 
        bwa index {params.args} {input} -p {params.index_path} 
        """

rule bwa_mem:
    input:
        reads = lambda wildcards: alignment_input(wildcards.sample),
        genomeIndex = multiext(f"{RESULTS}/bwa-index/{genome}.", "amb", "ann", "pac", "sa", "bwt")
    output:
        temp(RESULTS + "/bwa-mem/{sample}.bam")
    conda:
        "../envs/align.yml"
    params:
        genome = genome,
        args = config["bwa-mem"]["args"],
        extra = config["bwa-mem"]["extra"],
        index_path= f"{RESULTS}/bwa-index/{genome}",
    threads:
        int(config['bwa-mem']['threads'])
    log:
        LOGS + "/bwa-mem/{sample}.log"
    benchmark:
        BENCHMARKS + "/bwa-mem/{sample}.txt"
    resources:
        tmpdir=TEMP,
        mem_per_cpu=2096,
        cpus_per_task=int(config["bwa-mem"]["threads"])
    shell:
        """
        exec > {log} 2>&1
        bwa mem -t {threads} {params.args} {params.extra} {params.index_path} {input.reads} | samtools sort -@ {threads} -o {output} -
        """

rule buildBWA2Index:
    input:
        f"{RESOURCES}/genomes/{genome}.fa.gz"
    output:
        multiext(f"{RESULTS}/bwa2-index/{genome}.", "amb", "ann", "pac")
    conda:
        "../envs/align.yml"
    params:
        index_path = f"{RESULTS}/bwa2-index/{genome}",
        dir = f"{RESULTS}/bwa-index/",
        args = config["bwa-index"]["args"]
    log:
        f"{LOGS}/bwa2-index/{genome}.log"
    benchmark:
        f"{BENCHMARKS}/bwa2-index/{genome}.txt"
    resources:
        tmpdir=TEMP,
        mem_mb=20000,
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.dir} 
        bwa-mem2 index {params.args} -p {params.index_path} {input} 
        """

rule bwa_mem2:
    input:
        reads = lambda wildcards: alignment_input(wildcards.sample),
        genomeIndex = multiext(f"{RESULTS}/bwa2-index/{genome}.", "amb", "ann", "pac", "0123", "bwt.2bit.64")
    output:
        temp(RESULTS + "/bwa-mem2/{sample}.bam")
    conda:
        "../envs/align.yml"
    params:
        genome = genome,
        args = config["bwa-mem2"]["args"],
        extra = config["bwa-mem2"]["extra"],
        index_path= f"{RESULTS}/bwa2-index/{genome}",
    threads:
        int(config['bwa-mem2']['threads'])
    log:
        LOGS + "/bwa-mem2/{sample}.log"
    benchmark:
        BENCHMARKS + "/bwa-mem2/{sample}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_task=int(config['bwa-mem2']['threads']),
        mem_per_cpu=2096
    shell:
        """
        exec > {log} 2>&1
        bwa-mem2 mem -t {threads} {params.args} {params.extra} {params.index_path} {input.reads} | samtools sort -@ {threads} -o {output} -
        """

rule buildStarIndex:
    input:
        f"{RESOURCES}/genomes/{genome}.fa"
    output:
        multiext(RESULTS + "/star-index/SA", "", "index")
    params:
        genome_path = f"{RESOURCES}/genomes/{genome}.fa",
        result_path = f"{RESULTS}/star-index",
        args = config["STAR"]["args"]
    conda:
        "../envs/align.yml"
    threads:
        config["STAR"]["threads"]
    log:
        f"{LOGS}/star-index/{genome}.log"
    benchmark:
        f"{BENCHMARKS}/star-index/{genome}.txt"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        rm -rf {resources.tmpdir}/star-index
        STAR --outTmpDir {resources.tmpdir}/star-index --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.result_path} --genomeFastaFiles {params.genome_path} {params.args}
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
        STAR --outTmpDir "{resources.tmpdir}/STAR-${{uiidgen}}" --readFilesType Fastx --runThreadN {threads} --genomeDir {params.index_path} --readFilesIn {input.reads} {params.args} {params.extra} --outStd SAM | samtools sort -@ {threads} -o {output} -
        """
