rule unzip:
    input:
        "{file}.fa.gz"
    output:
        "{file}.fa"
    wildcard_constraints:
        file = r"^([\/A-Za-z0-9])*"
    shell:
        "gzip -dk {input} -f"


rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    threads:
        2
    conda:
        "../envs/utils.yml"
    log:
        "logs/samtools-index/{sample}.log"
    shell:
        """
        exec > {log} 2>&1
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

rule picard_MarkDuplicates:
    input:
        aligned = f"results/{config['aligner']}/{{sample}}.bam",
        aligned_index= "results/picard-MarkDuplicates/{sample}.bam.bai",
    output:
        sorted = temp("temp/{sample}_sorted.bam"),
        marked = multiext("results/picard-MarkDuplicates/{sample}", ".metrics.txt", ".bam"),
    params:
        paired_end = config['paired_end'],
        genome = config['genome'],
        aligner = config['aligner'],
        args = config['MarkDuplicates']['args']
    conda:
        "../envs/utils.yml"
    log:
        "logs/picard-MarkDuplicates/{sample}.log"
    shell:
        """
        exec > {log} 2>&1
        picard SortSam -I {input.aligned}  -O {output.sorted} --SO coordinate
        picard MarkDuplicates -I {output.sorted} -O results/picard-MarkDuplicates/{wildcards.sample}.bam -M results/picard-MarkDuplicates/{wildcards.sample}.metrics.txt {params.args}
        """

rule samtools_markdup:
    input:
        f"results/{config['aligner']}/{{sample}}.bam"
    output:
        "results/samtools-markdup/{sample}.bam"
    conda:
        "../envs/utils.yml"
    params:
        args = config['markdup']['args']
    threads:
        8
    log:
        "logs/samtools-markdup/{sample}.log"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p results/samtools-markdup
        samtools collate -@ {threads} -O -u {input} | samtools fixmate -@ {threads} -m -u - - | samtools sort -@ {threads} -u - | samtools markdup {params.args} -@ {threads} - {output}
        """