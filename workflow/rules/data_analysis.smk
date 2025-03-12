rule deeptools_bamCoverage:
    input:
        bam = RESULTS + "/" + config['duplicate_processor'] + "/{sample}.bam",
        bam_index = RESULTS + "/" + config['duplicate_processor'] + "/{sample}.bam.bai",
    output:
        RESULTS + "/deeptools-bamCoverage/{sample}.bw"
    conda:
        "../envs/data_analysis.yml"
    log:
        LOGS + "/bamCoverage/{sample}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards,threads: threads
    threads:
        12
    shell:
        """
        exec > {log} 2>&1
        bamCoverage -p {threads} -b {input.bam} -o {output}
        """

rule deeptools_computeMatrix:
    input:
        beds = lambda wildcards:f"{RESULTS}/{config['peak_caller']}/{wildcards.sample}_peaks.bed",
        bigwigs = lambda wildcards: f"{RESULTS}/deeptools-bamCoverage/{wildcards.sample}.bw"
    output:
        RESULTS + "/deeptools/{sample}_matrix.gz"
    wildcard_constraints:
        replicate = r"[0-9]"
    conda:
        "../envs/data_analysis.yml"
    params:
        mode = "reference-point",
        args = config["computeMatrix"]["args"],
        outdir = f"{RESULTS}/deeptools"
    threads:
        int(config["computeMatrix"]["threads"])
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards,threads: threads
    shell:
        """
        mkdir -p {params.outdir}
        computeMatrix {params.mode} -p {threads} -S {input.bigwigs} -R {input.beds} -o {output} {params.args}
        """

rule deeptools_plotHeatMap:
    input:
        RESULTS + "/deeptools/{sample}_matrix.gz"
    output:
        RESULTS + "/deeptools/{sample}_heatmap.png"
    conda:
        "../envs/data_analysis.yml"
    resources:
        tmpdir=TEMP
    shell:
        """
        plotHeatmap -m {input} -o {output}
        """
rule deeptools_plotProfile:
    input:
        RESULTS + "/deeptools/{sample}_matrix.gz"
    output:
        RESULTS + "/deeptools/{sample}_profile.png"
    conda:
        "../envs/data_analysis.yml"
    resources:
        tmpdir=TEMP
    shell:
        """
        plotProfile -m {input} -o {output}
        """

rule plot_genome_track:
    input:
        bed = lambda wildcards: [f"{RESULTS}/{config['peak_caller']}/{wildcards.sample}_peaks.bed"],
        bigwig = lambda wildcards: f"{RESULTS}/deeptools-bamCoverage/{wildcards.sample}.bw"
    output:
        tracks = temp(RESULTS + "/pyGenomeTracks/{sample}_tracks.ini"),
        plot = RESULTS + "/pyGenomeTracks/{sample}.png"
    params:
        region = config["plot_regions"],
        args = config["pyGenomeTracks"]["args"],
        peak_type = lambda wildcards: sfs.get_sample_entry_by_file_name(wildcards.sample)["peak_type"],
        bigwig_options = config["pyGenomeTracks"]["bigwig_options"],
        bed_options = config["pyGenomeTracks"]["bed_options"],
    conda:
        "../envs/data_analysis.yml"
    log:
        LOGS + "/pyGenomeTracks/{sample}.log"
    resources:
        tmpdir=TEMP
    script:
        '../scripts/pyGenomeTracks.py'

rule bedtools_intersect:
    input:
        a = lambda wildcards: get_consensus_peak_input(wildcards.sample)[0],
        b = lambda wildcards: get_consensus_peak_input(wildcards.sample)[1:]
    output:
        RESULTS + "/bedtools/{sample}.consensusPeak"
    conda:
        "../envs/data_analysis.yml"
    log:
        LOGS + "/bedtools-intersect/{sample}.log"
    benchmark:
        BENCHMARKS + "/bedtools-intersect/{sample}.txt"
    resources:
        tmpdir=TEMP,
    shell:
        '''
        exec > {log} 2>&1
        bedtools intersect -a {input.a} -b {input.b} -wa > {output}
        '''

rule homer_setup:
    output:
        touch(f"{RESULTS}/homer/homer_setup.done")
    conda:
        "../envs/data_analysis.yml"
    params:
        genome=genome
    resources:
        tmpdir=TEMP,
        cpus_per_task=lambda wildcards, threads: threads
    shell:
        '''
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/configureHomer.pl -install {params.genome} 
        '''


rule homer_annotate_peaks:
    input:
        lambda wildcards: f"{RESULTS}/bedtools/{wildcards.sample}.consensusPeak",
        f"{RESULTS}/homer/homer_setup.done"
    output:
        RESULTS + "/homer/{sample}_annotate.txt"
    conda:
        "../envs/data_analysis.yml"
    params:
        outdir = RESULTS + "/homer",
        genome = genome
    log:
        LOGS + "/homer/{sample}_annotate.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task = lambda wildcards, threads: threads
    shell:
        '''
        exec > {log} 2>&1
        mkdir -p {params.outdir} 
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/configureHomer.pl -install {params.genome} 
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/annotatePeaks.pl {input} {params.genome} > {output}
        '''

rule homer_find_motifs_genome:
    input:
        RESULTS + "/homer/{sample}_annotate.txt"
    output:
        multiext(RESULTS + "/homer/{sample}/", "homerResults.html", "knownResults.html")
    conda:
        "../envs/data_analysis.yml"
    params:
        outdir = RESULTS + "/homer",
        genome = config['genome'],
        args = config["homer"]["args"]
    threads:
        int(config["homer"]["threads"])
    log:
        LOGS + "/homer/{sample}_findMotifs.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards,threads: threads
    shell:
        '''
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/findMotifsGenome.pl {input} {params.genome} {params.outdir}/{wildcards.sample} -p {threads} {params.args}
        '''

rule plot_annotated_peaks:
    input:
        RESULTS + "/homer/{sample}_annotate.txt"
    output:
        distribution_plot = RESULTS + "/plots/{sample}_distribution.png",
        gene_distribution_plot =  RESULTS + "/plots/{sample}_genes.png"
    params:
        threshold_fraction = 0.01
    resources:
        mem_mb=200,
    conda:
        "../envs/data_analysis.yml"
    script:
        "../scripts/plot_annotated_peaks.py"
