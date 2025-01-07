import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import SampleFileScripts
sfs = SampleFileScripts(config)

def get_trim_input(sample: str, path: str, ext: str) -> list[str]:
    if sfs.is_paired_end():
        input = [f"{path}/{sample}_1.{ext}", f"{path}/{sample}_2.{ext}"]
    else:
        input = [f"{path}/{sample}.{ext}"]
    return input

rule trim_galore:
    input:
        lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq")
    output:
        expand("results/trim_galore/{sample}{read}.fastq", read=["_1", "_2"] if sfs.is_paired_end() else [""], allow_missing=True),
    conda:
        "../envs/trim.yml"
    params:
        args = config["trim_galore"]["args"],
        output_dir = "results/trim_galore",
        paired_end = config["paired_end"]
    threads:
        int(config["trim_galore"]["threads"])
    log:
        "logs/trim_galore/{sample}"
    shell:
        """
        exec > {log} 2>&1
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            read -a output_array <<< "{output}"
            out1="${{output_array[0]}}"
            out2="${{output_array[1]}}"
            trim_galore --paired -j {threads} -o {params.output_dir} --basename {wildcards.sample} {params.args} {input}
            mv {params.output_dir}/{wildcards.sample}_val_1.fq $out1
            mv {params.output_dir}/{wildcards.sample}_val_2.fq $out2
        else
            trim_galore -j {threads} -o {params.output_dir} --basename {wildcards.sample} {params.args} {input}
            mv {params.output_dir}/{wildcards.sample}_trimmed.fq {output} 
        fi
        """

rule cutadapt:
    input:
        lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq")
    output:
        expand("results/cutadapt/{sample}{read}.fastq", read=["_1", "_2"] if sfs.is_paired_end() else[""], allow_missing=True)
    conda:
        "../envs/cutadapt.yml"
    params:
        args = config["cutadapt"]["args"],
        paired_end = config["paired_end"]
    threads:
        int(config["cutadapt"]["threads"])
    log:
        "logs/cutadapt/{sample}.log"
    shell:
        '''
        exec > {log} 2>&1
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            read -a output_array <<< "{output}"
            out1="${{output_array[0]}}"
            out2="${{output_array[1]}}"
            cutadapt -o $out1 -p $out2 {params.args} {input}
        else
            cutadapt -j {threads} -o {output} {params.args} {input}
        fi
        '''

rule fastp:
    input:
        samples = lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq")
    output:
        expand("results/fastp/{sample}{read}.fastq", read=["_1", "_2"] if sfs.is_paired_end() else[""], allow_missing=True)
    conda:
        "../envs/trim.yml"
    params:
        args = config["fastp"]["args"],
        paired_end = config["paired_end"]
    threads:
        8
    log:
        "logs/fastp/{sample}.log"
    shell:
        '''
        exec > {log} 2>&1
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            read -a output_array <<< "{output}"
            out1="${{output_array[0]}}"
            out2="${{output_array[1]}}"
            fastp -j results/fastp/{wildcards.sample}.json -h results/fastp/{wildcards.sample}.html -i {input.samples[0]} -I {input.samples[1]} -o $out1 -O $out2 {params.args}
        else
            fastp -w {threads} -j results/fastp/{wildcards.sample}.json -h results/fastp/{wildcards.sample}.html -i {input} -o {output} {params.args}
        fi
        '''
