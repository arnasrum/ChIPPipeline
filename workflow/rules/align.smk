rule build_bowtie2_index:
    input:
        f"{RESOURCES}/genomes/{{genome}}.fa.gz"
    output:
        multiext(f"{RESULTS}/bowtie2-build/{{genome}}.", "1.bt2", "2.bt2", "3.bt2", "4.bt2"),
    conda:
        "../envs/align.yml"
    params:
        args = config["bowtie2-build"]["args"],
        path = RESULTS + "/bowtie2-build"
    threads:
        int(config["indexing_threads"])
    log:
        f"{LOGS}/bowtie2-build/{{genome}}.log"
    benchmark:
        f"{BENCHMARKS}/bowtie2-build/{{genome}}.txt"
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads,
    shell:
        '''
        exec > {log} 2>&1
        bowtie2-build --threads {threads} {params.args} {input} {params.path}/{wildcards.genome}
        '''

rule bowtie2:
    input:
        genome_index = lambda wildcards: multiext(f"{RESULTS}/bowtie2-build/{sfs.get_sample_genome(wildcards.sample)}.", "1.bt2", "2.bt2", "3.bt2", "4.bt2"),
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        RESULTS + "/bowtie2/{sample}.bam"
    conda:
        "../envs/align.yml"
    params:
        args = config["bowtie2"]["args"],
        genome = lambda wildcards: sfs.get_sample_genome(wildcards.sample),
        paired_end = config["paired_end"],
        index_path= lambda wildcards,input: os.path.commonpath(input.genome_index)
    threads:
        int(config["aligning_threads"])
    log:
        LOGS + "/bowtie2/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/bowtie2/{sample}.txt", int(config["benchmark_repeat_align"]))
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards, threads: threads
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
        bowtie2 -q --threads {threads} -x {params.index_path} $inputOptions {params.args} | samtools sort -@ {threads} -o {output} -
        '''
rule build_bwa_index:
    input:
        f"{RESOURCES}/genomes/{{genome}}.fa.gz"
    output: 
        multiext(f"{RESULTS}/bwa-index/{{genome}}.", "amb", "ann", "pac", "sa", "bwt")
    conda:
        "../envs/align.yml"
    params:
        index_path = f"{RESULTS}/bwa-index/{{genome}}",
        dir = f"{RESULTS}/bwa-index/",
        args = config["bwa-index"]["args"]
    log:
        f"{LOGS}/bwa-index/{{genome}}.log"
    benchmark:
        f"{BENCHMARKS}/bwa-index/{{genome}}.txt"
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
        genome_index = lambda wildcards: multiext(f"{RESULTS}/bwa-index/{sfs.get_sample_genome(wildcards.sample)}.", "amb", "ann", "pac", "sa", "bwt")
    output:
        temp(RESULTS + "/bwa-mem/{sample}.bam")
    conda:
        "../envs/align.yml"
    params:
        genome = lambda wildcards: sfs.get_sample_genome(wildcards.sample),
        args = config["bwa-mem"]["args"],
        index_path= lambda wildcards, input: f"{os.path.commonpath(input.genome_index)}/{sfs.get_sample_genome(wildcards.sample)}",
    threads:
        int(config['aligning_threads'])
    log:
        LOGS + "/bwa-mem/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/bwa-mem/{sample}.txt", int(config["benchmark_repeat_align"]))
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards, threads: threads
    shell:
        """
        exec > {log} 2>&1
        bwa mem -t {threads} {params.args} {params.index_path} {input.reads} | samtools sort -@ {threads} -o {output} -
        """

rule build_bwa2_index:
    input:
        f"{RESOURCES}/genomes/{{genome}}.fa.gz"
    output:
        multiext(f"{RESULTS}/bwa2-index/{{genome}}.", "amb", "ann", "pac", "bwt.2bit.64", "0123")
    conda:
        "../envs/align.yml"
    params:
        index_path = f"{RESULTS}/bwa2-index/{{genome}}",
        dir = f"{RESULTS}/bwa-index/",
        args = config["bwa-index"]["args"]
    log:
        f"{LOGS}/bwa2-index/{{genome}}.log"
    benchmark:
        f"{BENCHMARKS}/bwa2-index/{{genome}}.txt"
    resources:
        tmpdir=TEMP,
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.dir} 
        bwa-mem2 index {params.args} -p {params.index_path} {input} 
        """

rule bwa_mem2:
    input:
        reads = lambda wildcards: alignment_input(wildcards.sample),
        genome_index = lambda wildcards: multiext(f"{RESULTS}/bwa2-index/{sfs.get_sample_genome(wildcards.sample)}.", "amb", "ann", "pac", "bwt.2bit.64", "0123")
    output:
        temp(RESULTS + "/bwa-mem2/{sample}.bam")
    conda:
        "../envs/align.yml"
    params:
        genome = lambda wildcards: sfs.get_sample_genome(wildcards.sample),
        args = config["bwa-mem2"]["args"],
        index_path = lambda wildcards, input: os.path.commonpath(input.genome_index)
    threads:
        int(config['aligning_threads'])
    log:
        LOGS + "/bwa-mem2/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/bwa-mem2/{sample}.txt", int(config["benchmark_repeat_align"]))
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards, threads: threads
    shell:
        """
        exec > {log} 2>&1
        bwa-mem2 mem -t {threads} {params.args} {params.index_path} {input.reads} | samtools sort -@ {threads} -o {output} -
        """

rule build_STAR_index:
    input:
        f"{RESOURCES}/genomes/{{genome}}.fa"
    output:
        multiext(RESULTS + "/star_index/{genome}/SA", "", "index")
    params:
        genome_path = f"{RESOURCES}/genomes/{{genome}}.fa",
        result_path = lambda wildcards: f"{RESULTS}/star_index/{wildcards.genome}",
    conda:
        "../envs/align.yml"
    threads:
        config["indexing_threads"]
    log:
        f"{LOGS}/star_index/{{genome}}.log"
    benchmark:
        f"{BENCHMARKS}/star_index/{{genome}}.txt"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        rm -rf {resources.tmpdir}/star-index
        STAR --outTmpDir {resources.tmpdir}/star_index --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.result_path} --genomeFastaFiles {params.genome_path}
        """ 

rule STAR:
    input:
        genome_index = lambda wildcards: multiext(RESULTS + f"/star_index/{sfs.get_sample_genome(wildcards.sample)}/SA", "", "index"),
        reads = lambda wildcards: alignment_input(wildcards.sample)
    output:
        temp(RESULTS + "/STAR/{sample}.bam")
    conda:
        "../envs/align.yml"
    threads:
        config["aligning_threads"]
    params:
        args = config["star"]["args"],
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