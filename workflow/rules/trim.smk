ruleorder: trim_galore_PE > trim_galore_SE

rule trim:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        [f"{RESULTS}/{config['trimmer']}/{{sample}}_1.fastq", f"{RESULTS}/{config['trimmer']}/{{sample}}_2.fastq"]
         if sfs.is_paired_end() else
        [f"{RESULTS}/{config['trimmer']}/{{sample}}.fastq"]
    conda:
        "../envs/trim.yml" if not config["trimmer"] == "cutadapt" else "../envs/cutadapt.yml"
    params:
        args=config[config['trimmer']]["args"],
    threads:
        int(config["trim_threads"])
    log:
        f"{LOGS}/{config['trimmer']}/{{sample}}.log"
    benchmark:
        repeat(f"{BENCHMARKS}/{config['trimmer']}/{{sample}}.txt",config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
    script:
        "../scripts/trimmers/trim.py"


rule trim_galore_PE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        out1 = RESULTS + "/trim_galore/{sample}_1.fastq",
        out2 = RESULTS + "/trim_galore/{sample}_2.fastq"
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
        repeat(BENCHMARKS + "/trim_galore/{sample}.txt", config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=int(config["trim_galore"]["threads"]),
    shell:
        """
        exec > {log} 2>&1
        trim_galore --paired -j {threads} -o {params.output_dir} --basename {wildcards.sample} {params.args} {input}
        mv {params.output_dir}/{wildcards.sample}_val_1.fq {output.out1} 
        mv {params.output_dir}/{wildcards.sample}_val_2.fq {output.out2} 
        """

rule trim_galore_SE:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        RESULTS + "/trim_galore/{sample}.fastq"
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
        repeat(BENCHMARKS + "/trim_galore/{sample}.txt", config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=int(config["trim_galore"]["threads"]),
    shell:
        """
        exec > {log} 2>&1
        trim_galore -j {threads} -o {params.output_dir} --basename {wildcards.sample} {params.args} {input}
        mv {params.output_dir}/{wildcards.sample}_trimmed.fq {output} 
        """

rule cutadapt:
    input:
        lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/cutadapt/{sample}{extension}", extension=fastq_file_extensions, allow_missing=True))
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
        repeat(BENCHMARKS + "/cutadapt/{sample}.txt", config["benchmark_repeat_trim"])
    resources:
        tmpdir=TEMP,
        cpus_per_task=int(config["cutadapt"]["threads"]),
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