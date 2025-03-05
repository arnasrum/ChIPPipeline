

ruleorder: trim_galore_PE > trim_galore_SE

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

rule fastp:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample)
    output:
        temp(expand(RESULTS + "/fastp/{sample}{extension}", extension=fastq_file_extensions, allow_missing=True))
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
        cpus_per_task=int(config["fastp"]["threads"]),
    benchmark:
        repeat(BENCHMARKS + "/fastp/{sample}.txt", config["benchmark_repeat_trim"])
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

rule trimmomatic:
    input:
        samples = lambda wildcards: trimmer_input(wildcards.sample)
    output:
        samples = temp(expand(RESULTS + "/trimmomatic/{sample}{extension}", extension=fastq_file_extensions, allow_missing=True)),
    conda:
        "../envs/trim.yml"
    params:
        args = config["trimmomatic"]["args"],
        paired_end = config["paired_end"],
        run_options = config["trimmomatic"]["run_options"]
    threads:
        int(config["trimmomatic"]["threads"])
    log:
        log = LOGS + "/trimmomatic/{sample}.log",
        summary= LOGS + "/trimmomatic/{sample}_summary.txt",
    resources:
        tmpdir=TEMP,
    benchmark:
        repeat(BENCHMARKS + "/trimmomatic/{sample}.txt", config["benchmark_repeat_trim"])
    script:
        "../scripts/trimmomatic.py"