rule unzip_genome:
    input:
        f"{RESOURCES}/genomes/{{genome}}.gz"
    output:
        temp(f"{RESOURCES}/genomes/{{genome}}")
    wildcard_constraints:
        genome = r"[.]*.(fastq|fa)"
    resources:
        tmpdir=TEMP
    shell:
        "gzip -dkf {input}"

rule unzip_sample:
    input:
        f"{RESULTS}/{config['trimmer']}/{{sample}}.gz"
    output:
        temp(f"{RESULTS}/{config['trimmer']}/{{sample}}")
    wildcard_constraints:
        sample = r"[.]*.fastq"
    resources:
        tmpdir=TEMP
    shell:
        "gzip -dkf {input}"

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

