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
        int(config['bowtie2-build']['threads'])
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
    threads:
        int(config['bowtie2']['threads'])
    log:
        LOGS + "/bowtie2/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/bowtie2/{sample}.txt", int(config["benchmark_repeat_align"]))
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, attempt: config['bowtie2']['mem_mb'] * attempt,
        runtime=lambda wildcards, attempt: config['bowtie2']['runtime'] * attempt
    script:
        "../scripts/tools/bowtie2.py"

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
        args = config["bwa_mem_index"]["args"]
    log:
        f"{LOGS}/bwa-index/{{genome}}.log"
    benchmark:
        f"{BENCHMARKS}/bwa-index/{{genome}}.txt"
    resources:
        tmpdir=TEMP,
        mem_mb= lambda wildcards,attempt: config['bwa_mem_index']['mem_mb'] * attempt,
        runtime= lambda wildcards,attempt: config['bwa_mem_index']['runtime'] * attempt
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
        temp(RESULTS + "/bwa_mem/{sample}.bam")
    conda:
        "../envs/align.yml"
    params:
        args = config["bwa_mem"]["args"],
    threads:
        int(config['bwa_mem']['threads'])
    log:
        LOGS + "/bwa_mem/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/bwa_mem/{sample}.txt", int(config["benchmark_repeat_align"]))
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards, threads: threads,
        mem_mb= lambda wildcards,attempt: config['bwa_mem']['mem_mb'] * attempt,
        runtime= lambda wildcards,attempt: config['bwa_mem']['runtime'] * attempt
    script:
        "../scripts/tools/bwa-mem.py"

rule build_bwa2_index:
    input:
        f"{RESOURCES}/genomes/{{genome}}.fa.gz"
    output:
        multiext(f"{RESULTS}/bwa2_index/{{genome}}.", "amb", "ann", "pac", "bwt.2bit.64", "0123")
    conda:
        "../envs/align.yml"
    params:
        index_path = f"{RESULTS}/bwa2_index/{{genome}}",
        dir = f"{RESULTS}/bwa2_index/",
        args = config["bwa_mem2_index"]["args"]
    log:
        f"{LOGS}/bwa2_index/{{genome}}.log"
    benchmark:
        f"{BENCHMARKS}/bwa2_index/{{genome}}.txt"
    resources:
        tmpdir=TEMP,
        mem_mb= lambda wildcards,attempt: config['bwa_mem2_index']['mem_mb'] * attempt,
        runtime= lambda wildcards,attempt: config['bwa_mem2_index']['runtime'] * attempt
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.dir} 
        bwa-mem2 index {params.args} -p {params.index_path} {input} 
        """

rule bwa_mem2:
    input:
        reads = lambda wildcards: alignment_input(wildcards.sample),
        genome_index = lambda wildcards: multiext(f"{RESULTS}/bwa2_index/{sfs.get_sample_genome(wildcards.sample)}.", "amb", "ann", "pac", "bwt.2bit.64", "0123")
    output:
        temp(RESULTS + "/bwa_mem2/{sample}.bam")
    conda:
        "../envs/align.yml"
    params:
        args = config["bwa_mem2"]["args"],
    threads:
        int(config['bwa_mem2']['threads'])
    log:
        LOGS + "/bwa_mem2/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/bwa_mem2/{sample}.txt", int(config["benchmark_repeat_align"]))
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards, threads: threads,
        mem_mb= lambda wildcards,attempt: config['bwa_mem2']['mem_mb'] * attempt,
        runtime= lambda wildcards,attempt: config['bwa_mem2']['runtime'] * attempt
    script:
        "../scripts/tools/bwa-mem2.py"

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
        int(config['star_index']['threads'])
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
        reads = lambda wildcards: [sample.replace(".gz", "") for sample in alignment_input(wildcards.sample)]
    output:
        temp(RESULTS + "/STAR/{sample}.bam")
    conda:
        "../envs/align.yml"
    threads:
        int(config['star']['threads'])
    params:
        args = config["star"]["args"],
        index_path = f"{RESULTS}/star-index",
        output_path = f"{RESULTS}/STAR",
        read_group= '@RG\tID:{sample}\tSM:{sample}'
    log:
        LOGS + "/STAR/{sample}.log"
    benchmark:
        repeat(BENCHMARKS + "/STAR/{sample}.txt", config["benchmark_repeat_align"])
    resources:
        tmpdir=TEMP
    script:
        "../scripts/tools/star.py"
