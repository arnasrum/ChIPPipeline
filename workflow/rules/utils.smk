
ruleorder: PicardBuildBamIndex > samtools_index

rule unzip_genome:
    input:
        RESOURCES + "/genomes/{genome}.gz"
    output:
        RESOURCES + "/genomes/{genome}"
    wildcard_constraints:
        genome = r"^(.*).(fa|fasta)$"
    resources:
        tmpdir=TEMP
    params:
        args = config["gzip"]["args"]
    shell:
        "gzip {params.args} -dkf {input}"

rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
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
        samtools sort -@ {threads} {input} | samtools index -@ {threads} - -o {output}
        """

rule PicardBuildBamIndex:
    input:
        "{sample}.bam"
    output:
        sorted = temp(f"{TEMP}/{{sample}}.bam"),
        index = "{sample}.bam.bai"
    log:
        LOGS + "/BuildBamIndex/{sample}.log"
    conda:
        "../envs/utils.yml"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        picard SortSam -I {input} -O {output.sorted} -SO coordinate --TMP_DIR {resources.tmpdir}
        picard BuildBamIndex -I {output.sorted} -O {output.index}
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
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output}"

