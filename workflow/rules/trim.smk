OUTPUTDIRS = {"trimgalore": "results/trim_galore", "cutadapt": "results/cutadapt", "fastp": "results/fastp"}

rule trimgalore_pe:
    input:
        "resources/reads/{id}_1.fastq",  
        "resources/reads/{id}_2.fastq"
    output:
        out1 = f"{OUTPUTDIRS["trimgalore"]}/{{id}}_1.fastq",
        out2 = f"{OUTPUTDIRS["trimgalore"]}/{{id}}_2.fastq"
    conda:
        "../envs/trim.yml"
    params:
        name = f"{{id}}",
        args = config["trimgalore"]["args"],
        outputdir = OUTPUTDIRS["trimgalore"]

    shell:
        """
        trim_galore --paired --no_report_file -o {params.outputdir} --basename {params.name} {params.args} {input}
        mv {params.outputdir}/{wildcards.id}_val_1.fq {output.out1} 
        mv {params.outputdir}/{wildcards.id}_val_2.fq {output.out2}
        """

rule cutadapt_pe:
    input:
        "resources/reads/{id}_1.fastq",
        "resources/reads/{id}_2.fastq"
    output:
        out1 = f"{OUTPUTDIRS["cutadapt"]}/{{id}}_1.fastq",
        out2 = f"{OUTPUTDIRS["cutadapt"]}/{{id}}_2.fastq"
    conda:
        "../envs/cutadapt.yml"
    params:
        args = config["cutadapt"]["args"]
    shell:
        '''
        cutadapt -o {output.out1} -p {output.out2} {params.args} {input}
        '''

rule fastp_pe:
    input:
        read1 = "resources/reads/{id}_1.fastq",
        read2 = "resources/reads/{id}_2.fastq"
    output:
        out1 = f"{OUTPUTDIRS["fastp"]}/{{id}}_1.fastq",
        out2 = f"{OUTPUTDIRS["fastp"]}/{{id}}_2.fastq"
    conda:
        "../envs/trim.yml"
    params:
        args = config["fastp"]["args"]
    shell:
        '''
        fastp -j results/fastp/fastp.json -h results/fastp/fastp.html -i {input.read1} -I {input.read2} -o {output.out1} -O {output.out2} {params.args}
        '''
