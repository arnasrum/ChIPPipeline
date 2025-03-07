
rule build_index:
    input:
        f"{RESOURCES}/genomes/{genome}.fa.gz"
    output:
        index_files
    conda:
        "../envs/align.yml"
    params:
        args=config["index_args"],
    threads:
        int(config["indexing_threads"])
    log:
        f"{LOGS}/{aligner.get_name()}_index/{genome}.log"
    benchmark:
        f"{BENCHMARKS}/{aligner.get_name()}_index/{genome}.txt"
    resources:
        tmpdir=TEMP,
    script:
        "../scripts/aligners/build_index.py"

rule align:
    input:
        index = index_files,
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        f"{RESULTS}/{aligner_name}/{{sample}}.bam"
    conda:
        "../envs/align.yml"
    params:
        args=config[aligner_name]["args"],
    threads:
        int(config["aligning_threads"])
    log:
        f"{LOGS}/{aligner_name}/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/{aligner_name}/{{sample}}.txt"
    resources:
        tmpdir=TEMP,
    script:
        "../scripts/aligners/align.py"

rule build_bwa_index:
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
        index_path= f"{RESULTS}/bwa-index/{genome}",
    log:
        LOGS + "/bwa-mem/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/bwa-mem/{sample}.txt", int(config["benchmark_repeat_align"]))
    resources:
        tmpdir=TEMP,
    shell:
        """
        exec > {log} 2>&1
        bwa mem -t {threads} {params.args} {params.index_path} {input.reads} | samtools sort -@ {threads} -o {output} -
        """

rule build_STAR_index:
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
    #threads:
        #config["STAR"]["threads"]
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
    #threads:
        #config["STAR"]["threads"]
    params:
        args = config["STAR"]["args"],
        index_path = f"{RESULTS}/star-index",
        output_path = f"{RESULTS}/STAR"
    log:
        LOGS + "/STAR/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/STAR/{sample}.txt", config["benchmark_repeat_align"])
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        uuid=$(python3 -c "import uuid; print(uuid.uuid4())")
        STAR --outTmpDir "{resources.tmpdir}/STAR-${{uuid}}" --outFileNamePrefix {params.output_path}/{wildcards.sample}_ --readFilesType Fastx --runThreadN {threads} --genomeDir {params.index_path} --readFilesIn {input.reads} {params.args} --outStd SAM | samtools sort -@ {threads} -o {output} -
        """
