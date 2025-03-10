
ruleorder: PicardBuildBamIndex > samtools_index

rule unzip:
    input:
        "{file}.gz"
    output:
        "{file}"
    wildcard_constraints:
        file = r"^(.*).(fastq|fq|fa|fasta)$"
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
        samtools index -@ {threads} - -o {output}
        """

rule PicardBuildBamIndex:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
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
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output}"

