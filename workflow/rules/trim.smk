import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import SampleFileScripts
sfs = SampleFileScripts(config)

RESULTS: str = config['results_path']
RESOURCES: str = config['resources_path']
LOGS: str = config['logs_path']
TEMP: str = config['temp_path']
BENCHMARKS: str = config['benchmarks_path']

def get_trim_input(sample: str) -> list[str]:
    if sfs.is_paired_end():
        input = [RESOURCES + f"/reads/{sample}_1.fastq", RESOURCES + f"/reads/{sample}_2.fastq"]
    else:
        input = [RESOURCES + f"/reads/{sample}.fastq"]
    return input

read_extention = ["_1", "_2"] if sfs.is_paired_end() else [""]
rule trim_galore:
    input:
        lambda wildcards: get_trim_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/trim_galore/{sample}{read}.fastq", read=read_extention, allow_missing=True))
    conda:
        "../envs/trim.yml"
    params:
        args = config["trim_galore"]["args"],
        output_dir = RESULTS + "/trim_galore",
        paired_end = config["paired_end"]
    threads:
        int(config["trim_galore"]["threads"])
    log:
        LOGS + "/trim_galore/{sample}.log"
    benchmark:
        BENCHMARKS + "/trim_galore/{sample}.txt"
    resources:
        tmpdir=TEMP,
        cpu_per_task=int(config["trim_galore"]["threads"]),
        mem_per_cpu=2096
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
        lambda wildcards: get_trim_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/cutadapt/{sample}{read}.fastq", read=read_extention, allow_missing=True))
    conda:
        "../envs/cutadapt.yml"
    params:
        args = config["cutadapt"]["args"],
        paired_end = config["paired_end"]
    threads:
        int(config["cutadapt"]["threads"])
    log:
        LOGS + "/cutadapt/{sample}.log"
    benchmark:
        BENCHMARKS + "/cutadapt/{sample}.txt"
    resources:
        tmpdir=TEMP,
        cpu_per_task=int(config["cutadapt"]["threads"]),
        mem_per_cpu=2096
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
        samples = lambda wildcards: get_trim_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/fastp/{sample}{read}.fastq", read=read_extention, allow_missing=True))
    conda:
        "../envs/trim.yml"
    params:
        args = config["fastp"]["args"],
        paired_end = config["paired_end"],
        path = RESULTS + "/fastp"
    threads:
        int(config["fastp"]["threads"])
    log:
        LOGS + "/fastp/{sample}.log"
    resources:
        tmpdir=TEMP,
        cpu_per_task=int(config["fastp"]["threads"]),
        mem_per_cpu=2096
    shell:
        '''
        exec > {log} 2>&1
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            read -a output_array <<< "{output}"
            out1="${{output_array[0]}}"
            out2="${{output_array[1]}}"
            fastp -j {params.path}/{wildcards.sample}.json -h {params.path}/{wildcards.sample}.html -i {input.samples[0]} -I {input.samples[1]} -o $out1 -O $out2 {params.args}
        else
            fastp -w {threads} -j {params.path}/{wildcards.sample}.json -h {params.path}/{wildcards.sample}.html -i {input} -o {output} {params.args}
        fi
        '''
