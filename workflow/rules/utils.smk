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

