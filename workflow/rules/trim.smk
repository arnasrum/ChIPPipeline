import sys
sys.path.append("workflow/scripts")
from getFileNames import getFileNames 

OUTPUTDIRS = {"trimgalore": "results/trim_galore", "cutadapt": "results/cutadapt", "fastp": "results/fastp"}

rule trimgalore:
    input:
        lambda wildcards: getFileNames(wildcards.sample, "resources/reads", "fastq", config)
    output:
        out1 = f"{OUTPUTDIRS["trimgalore"]}/{{sample}}_1.fastq",
        out2 = f"{OUTPUTDIRS["trimgalore"]}/{{sample}}_2.fastq" if config["paired_end"] else []
        #out = getFileNames("wildcards.sample", "results/trim_galore", "fastq", config)
    conda:
        "../envs/trim.yml"
    params:
        name = f"{{sample}}",
        args = config["trimgalore"]["args"],
        outputdir = OUTPUTDIRS["trimgalore"],
        paired_end = config["paired_end"]
    shell:
        """
        if [[ {params.paired_end} ]]; then
            trim_galore --paired --no_report_file -o {params.outputdir} --basename {params.name} {params.args} {input}
            mv {params.outputdir}/{wildcards.sample}_val_1.fq {output.out1} 
            mv {params.outputdir}/{wildcards.sample}_val_2.fq {output.out2}
        else
            trim_galore --no_report_file -o {params.outputdir} --basename {params.name} {params.args} {input}
            mv {params.outputdir}/{wildcards.sample}_val_1.fq {output.out1} 
        fi
        """

rule cutadapt:
    input:
        lambda wildcards: getFileNames(wildcards.sample, "resources/reads", "fastq", config)
    output:
        out1 = f"{OUTPUTDIRS["cutadapt"]}/{{sample}}_1.fastq",
        out2 = f"{OUTPUTDIRS["cutadapt"]}/{{sample}}_2.fastq" if config["paired_end"] else []
    conda:
        "../envs/cutadapt.yml"
    params:
        args = config["cutadapt"]["args"],
        paired_end = config["paired_end"]
    shell:
        '''
        if [[ {params.paired_end} ]]; then
            cutadapt -o {output.out1} -p {output.out2} {params.args} {input}
        else
            cutadapt -o {output.out1} {params.args} {input}
        fi
        '''

rule fastp:
    input:
        samples = lambda wildcards: getFileNames(wildcards.sample, "resources/reads", "fastq", config)
    output:
        out1 = f"{OUTPUTDIRS["fastp"]}/{{sample}}_1.fastq",
        out2 = f"{OUTPUTDIRS["fastp"]}/{{sample}}_2.fastq" if config["paired_end"] else []
    conda:
        "../envs/trim.yml"
    params:
        args = config["fastp"]["args"],
        paired_end = config["paired_end"],
    shell:
        '''
        if [[ {params.paired_end} ]]; then
            fastp -j results/fastp/fastp.json -h results/fastp/fastp.html -i {input.samples[0]} -I {input.samples[1]} -o {output.out1} -O {output.out2} {params.args}
        else
            fastp -j results/fastp/fastp.json -h results/fastp/fastp.html -i {input.samples} -o {output.out1} {params.args}
        fi
        '''
