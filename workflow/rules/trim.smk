import sys
sys.path.append("workflow/scripts")

def get_trim_input(sample: str, path: str, ext: str, config: dict) -> list[str]:
    samples = [f"{path}/{sample}_1.{ext}"]
    if config["paired_end"]: samples.append(f"resources/reads/{sample}_2.{ext}")
    return samples

OUTPUTDIRS = {"trimgalore": "results/trim_galore", "cutadapt": "results/cutadapt", "fastp": "results/fastp"}

rule trimgalore:
    input:
        lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq", config)
    output:
        out1 = f"{OUTPUTDIRS['trimgalore']}/{{sample}}_1.fastq",
        out2 = f"{OUTPUTDIRS['trimgalore']}/{{sample}}_2.fastq" if config["paired_end"] else []
    conda:
        "../envs/trim.yml"
    params:
        name = f"{{sample}}",
        args = config["trimgalore"]["args"],
        output_dir = OUTPUTDIRS["trimgalore"],
        paired_end = config["paired_end"]
    shell:
        """
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            trim_galore --paired --no_report_file -o {params.output_dir} --basename {params.name} {params.args} {input}
            mv {params.output_dir}/{wildcards.sample}_val_1.fq {output.out1} 
            mv {params.output_dir}/{wildcards.sample}_val_2.fq {output.out2}
        else
            trim_galore --no_report_file -o {params.output_dir} --basename {params.name} {params.args} {input}
            mv {params.output_dir}/{wildcards.sample}_val_1.fq {output.out1} 
        fi
        """

rule cutadapt:
    input:
        lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq", config)
    output:
        out1 = f"{OUTPUTDIRS['cutadapt']}/{{sample}}_1.fastq",
        out2 = f"{OUTPUTDIRS['cutadapt']}/{{sample}}_2.fastq" if config["paired_end"] else []
    conda:
        "../envs/cutadapt.yml"
    params:
        args = config["cutadapt"]["args"],
        paired_end = config["paired_end"]
    shell:
        '''
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            cutadapt -o {output.out1} -p {output.out2} {params.args} {input}
        else
            cutadapt -o {output.out1} {params.args} {input}
        fi
        '''

rule fastp:
    input:
        samples = lambda wildcards: get_trim_input(wildcards.sample,"resources/reads","fastq", config)
    output:
        out1 = f"{OUTPUTDIRS['fastp']}/{{sample}}_1.fastq",
        out2 = f"{OUTPUTDIRS['fastp']}/{{sample}}_2.fastq" if config["paired_end"] else []
    conda:
        "../envs/trim.yml"
    params:
        args = config["fastp"]["args"],
        paired_end = config["paired_end"]
    shell:
        '''
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            fastp -j results/fastp/fastp.json -h results/fastp/fastp.html -i {input.samples[0]} -I {input.samples[1]} -o {output.out1} -O {output.out2} {params.args}
        else
            fastp -j results/fastp/fastp.json -h results/fastp/fastp.html -i {input.samples} -o {output.out1} {params.args}
        fi
        '''
