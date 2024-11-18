
rule unzip:
    input:
        "{file}.fa.gz"
    output:
        "{file}.fa"
    wildcard_constraints:
        file = r"^([\/A-Za-z0-9])*"
    shell:
        "gzip -dk {input} -f"


rule indexBam:
    input:
        f"results/{config['aligner']}/{{sample}}.bam"
    output:
        "results/samtools-index/{sample}.bam.bai"
    threads:
        2
    conda:
        "../envs/utils.yml"
    shell:
        """
        mkdir -p results/samtools-index 
        samtools index -@ {threads} {input} -o {output}
        """

rule picardCreateGenomeSequenceDictionary:
    input:
        "resources/genomes/{genome}.fa"
    output:
        "resources/genomes/{genome}.dict"
    conda:
        "../envs/utils.yml"
    shell:
        "picard CreateSequenceDictionary -R resources/genomes/{wildcards.genome}.fa -O resources/genomes/{wildcards.genome}.dict"



rule picardMarkDuplicates:
    input:
        f"results/{config['aligner']}/{{sample}}.bam",
        f"resources/genomes/{config['genome']}.dict",
        f"resources/genomes/{config['genome']}.fa",
        "results/samtools-index/{sample}.bam.bai"
    output:
        temp("temp/{sample}_unmapped.sam"),
        temp("temp/{sample}_merged.bam"),
        "results/picard-MarkDuplicates/{sample}.metrics.txt",
        "results/picard-MarkDuplicates/{sample}.bam",
        "results/picard-MarkDuplicates/{sample}.bam.bai"
    params:
        paired_end = config['paired_end'],
        genome = config['genome'],
        aligner = config['aligner']
    conda:
        "../envs/utils.yml"
    shell:
        """
        inputOptions='-F1 resources/reads/{wildcards.sample}_1.fastq'
        if [[ {params.paired_end} ]]; then
            inputOptions+=' -F2 resources/reads/{wildcards.sample}_2.fastq'
        fi 
        mkdir -p temp 
        picard FastqToSam ${{inputOptions}} -O temp/{wildcards.sample}_unmapped.sam -SAMPLE_NAME {wildcards.sample} 
        picard MergeBamAlignment -ALIGNED results/{params.aligner}/{wildcards.sample}.bam -UNMAPPED temp/{wildcards.sample}_unmapped.sam -O temp/{wildcards.sample}_merged.bam -R resources/genomes/{params.genome}.fa
        mkdir -p results/picard-MarkDuplicates
        picard MarkDuplicates -I temp/{wildcards.sample}_merged.bam -O results/picard-MarkDuplicates/{wildcards.sample}.bam -M results/picard-MarkDuplicates/{wildcards.sample}.metrics.txt
        cp results/samtools-index/{wildcards.sample}.bam.bai results/picard-MarkDuplicates/{wildcards.sample}.bam.bai
        """
