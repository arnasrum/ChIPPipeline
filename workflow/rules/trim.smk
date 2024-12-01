import sys
sys.path.append("workflow/scripts")

def get_trim_input(sample: str, path: str, ext: str, config: dict) -> list[str]:
    if config["paired_end"]: return [f"{path}/{sample}_1.{ext}", f"{path}/{sample}_2.{ext}"]
    else: return [f"{path}/{sample}.{ext}"]

OUTPUTDIRS = {"trimgalore": "results/trim_galore", "cutadapt": "results/cutadapt", "fastp": "results/fastp"}

rule trimgalore:
    input:
        lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq", config)
    output:
        #out1 = f"{OUTPUTDIRS['trimgalore']}/{{sample}}_1.fastq",
        #out2 = f"{OUTPUTDIRS['trimgalore']}/{{sample}}_2.fastq" if config["paired_end"] else [],
        expand("results/trim_galore/{sample}{read}.fastq", read=["_1", "_2"] if config["paired_end"] else [""], allow_missing=True)
    conda:
        "../envs/trim.yml"
    params:
        name = f"{{sample}}",
        args = config["trimgalore"]["args"],
        output_dir = OUTPUTDIRS["trimgalore"],
        paired_end = config["paired_end"]
    threads:
        8
    shell:
        """
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            read -a output_array <<< "{output}"
            out1="${{output_array[0]}}"
            out2="${{output_array[1]}}"
            trim_galore --paired --no_report_file -j {threads} -o {params.output_dir} --basename {params.name} {params.args} {input}
            mv {params.output_dir}/{wildcards.sample}_val_1.fq $out1
            mv {params.output_dir}/{wildcards.sample}_val_2.fq $out2
        else
            trim_galore --no_report_file -j {threads} -o {params.output_dir} --basename {params.name} {params.args} {input}
            mv {params.output_dir}/{wildcards.sample}_trimmed.fq {output} 
        fi
        """

rule cutadapt:
    input:
        lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq", config)
    output:
        #out1 = f"{OUTPUTDIRS['cutadapt']}/{{sample}}_1.fastq",
        #out2 = f"{OUTPUTDIRS['cutadapt']}/{{sample}}_2.fastq" if config["paired_end"] else []
        expand("results/cutadapt/{sample}{read}.fastq", read=["_1", "_2"] if config["paired_end"] else[""], allow_missing=True)

    conda:
        "../envs/cutadapt.yml"
    params:
        args = config["cutadapt"]["args"],
        paired_end = config["paired_end"]
    shell:
        '''
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            read -a output_array <<< "{output}"
            out1="${{output_array[0]}}"
            out2="${{output_array[1]}}"
            cutadapt -o $out1 -p $out2 {params.args} {input}
        else
            cutadapt -o {output} {params.args} {input}
        fi
        '''

rule fastp:
    input:
        samples = lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq", config)
    output:
        #out1 = f"{OUTPUTDIRS['fastp']}/{{sample}}_1.fastq",
        #out2 = f"{OUTPUTDIRS['fastp']}/{{sample}}_2.fastq" if config["paired_end"] else []
        expand("results/fastp/{sample}{read}.fastq", read=["_1", "_2"] if config["paired_end"] else[""], allow_missing=True)
    conda:
        "../envs/trim.yml"
    params:
        args = config["fastp"]["args"],
        paired_end = config["paired_end"]
    shell:
        '''
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            read -a output_array <<< "{output}"
            out1="${{output_array[0]}}"
            out2="${{output_array[1]}}"
            fastp -j results/fastp/fastp.json -h results/fastp/fastp.html -i {input.samples[0]} -I {input.samples[1]} -o $out1 -O $out2 {params.args}
        else
            fastp -j results/fastp/fastp.json -h results/fastp/fastp.html -i {input.samples} -o {output} {params.args}
        fi
        '''
