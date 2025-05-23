rule unzip_genome:
    input:
        f"{RESOURCES}/genomes/{{genome}}.gz"
    output:
        temp(f"{RESOURCES}/genomes/{{genome}}")
    resources:
        tmpdir=TEMP
    shell:
        "gzip -dkf {input}"

rule unzip_sample:
    input:
        f"{RESULTS}/{pipeline_config.get_config_option('trimmer')}/{{sample}}.gz"
    output:
        temp(f"{RESULTS}/{pipeline_config.get_config_option('trimmer')}/{{sample}}")
    resources:
        tmpdir=TEMP
    shell:
        "gzip -dkf {input}"

rule picardCreateGenomeSequenceDictionary:
    input:
        f"{RESOURCES}/genomes/{{genome}}.fa"
    output:
        f"{RESOURCES}/genomes/{{genome}}.dict"
    conda:
        "../envs/utils.yml"
    resources:
        tmpdir=TEMP,
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output}"

