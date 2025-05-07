rule deeptools_bamCoverage:
    input:
        bam = f"{RESULTS}/{config['duplicate_processor']}/{{sample}}.bam",
        bam_index = f"{RESULTS}/{config['duplicate_processor']}/{{sample}}.bam.bai",
    output:
        f"{RESULTS}/deeptools-bamCoverage/{{sample}}.bw"
    conda:
        "../envs/data_analysis.yml"
    log:
        f"{LOGS}/bamCoverage/{{sample}}.log"
    threads:
        int(config['bamCoverage']['threads'])
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards,threads: threads,
        mem_mb= lambda wildcards,attempt: int(config['bamCoverage']['mem_mb']) * attempt,
        runtime= lambda wildcards,attempt: int(config['bamCoverage']['runtime']) * attempt,
    shell:
        """
        exec > {log} 2>&1
        bamCoverage -p {threads} -b {input.bam} -o {output}
        """

rule deeptools_computeMatrix:
    input:
        beds = lambda wildcards:f"{RESULTS}/{config['peak_caller']}/{wildcards.sample}_peaks.narrowPeak",
        bigwigs = lambda wildcards: f"{RESULTS}/deeptools-bamCoverage/{wildcards.sample}.bw"
    output:
        f"{RESULTS}/deeptools/{{sample}}_matrix.gz"
    wildcard_constraints:
        replicate = r"[0-9]"
    conda:
        "../envs/data_analysis.yml"
    params:
        mode = config["computeMatrix"]["mode"],
        args = config["computeMatrix"]["args"],
        outdir = f"{RESULTS}/deeptools"
    threads:
        int(config["computeMatrix"]["threads"])
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards, threads: threads,
        mem_mb= lambda wildcards, attempt: int(config['computeMatrix']['mem_mb']) * attempt,
        runtime= lambda wildcards, attempt: int(config['computeMatrix']['runtime']) * attempt,
    shell:
        """
        mkdir -p {params.outdir}
        computeMatrix {params.mode} {params.args} -p {threads} -S {input.bigwigs} -R {input.beds} -o {output}
        """

rule deeptools_plotHeatMap:
    input:
        f"{RESULTS}/deeptools/{{sample}}_matrix.gz"
    output:
        f"{RESULTS}/deeptools/{{sample}}_heatmap.png"
    conda:
        "../envs/data_analysis.yml"
    resources:
        tmpdir=TEMP,
        runtime= lambda wildcards,attempt: config['plotHeatMap']['runtime'] * attempt,
        mem_mb= lambda wildcards,attempt: config['plotHeatMap']['mem_mb'] * attempt,
    shell:
        """
        plotHeatmap -m {input} -o {output} --whatToShow 'heatmap and colorbar'
        """
rule deeptools_plotProfile:
    input:
        f"{RESULTS}/deeptools/{{sample}}_matrix.gz"
    output:
        f"{RESULTS}/deeptools/{{sample}}_profile.png"
    conda:
        "../envs/data_analysis.yml"
    resources:
        tmpdir=TEMP,
        runtime = lambda wildcards, attempt: config['plotProfile']['runtime'] * attempt,
        mem_mb = lambda  wildcards, attempt: config['plotProfile']['mem_mb'] * attempt,
    shell:
        """
        plotProfile -m {input} -o {output}
        """

rule pyGenomeTracks:
    input:
        bed = lambda wildcards: [f"{RESULTS}/{config['peak_caller']}/{wildcards.sample}_peaks.narrowPeak"],
        bigwig = lambda wildcards: f"{RESULTS}/deeptools-bamCoverage/{wildcards.sample}.bw"
    output:
        tracks = temp(f"{RESULTS}/pyGenomeTracks/{{sample}}_tracks.ini"),
        plot = f"{RESULTS}/pyGenomeTracks/{{sample}}.png"
    params:
        region = config["plot_region"],
        args = config["pyGenomeTracks"]["args"],
        peak_type = lambda wildcards: pipeline_config.get_sample_entry_by_file_name(wildcards.sample)["peak_type"],
        bigwig_options = config["pyGenomeTracks"]["bigwig_options"],
        bed_options = config["pyGenomeTracks"]["bed_options"],
    conda:
        "../envs/data_analysis.yml"
    log:
        f"{LOGS}/pyGenomeTracks/{{sample}}.log"
    resources:
        tmpdir=TEMP
    script:
        '../scripts/pyGenomeTracks.py'

rule bedtools_intersect:
    input:
        a = lambda wildcards: get_consensus_peak_input(wildcards.group)[0],
        b = lambda wildcards: get_consensus_peak_input(wildcards.group)[1:]
    output:
        f"{RESULTS}/bedtools/{{group}}.consensusPeak"
    conda:
        "../envs/data_analysis.yml"
    log:
        f"{LOGS}/bedtools-intersect/{{group}}.log"
    benchmark:
        f"{BENCHMARKS}/bedtools-intersect/{{group}}.txt"
    resources:
        tmpdir=TEMP,
    shell:
        '''
        exec > {log} 2>&1
        bedtools intersect -a {input.a} -b {input.b} -wa > {output}
        '''

rule homer_annotate_peaks:
    input:
        lambda wildcards: f"{RESULTS}/bedtools/{wildcards.group}.consensusPeak",
    output:
        f"{RESULTS}/homer/{{group}}_annotate.txt"
    conda:
        "../envs/data_analysis.yml"
    params:
        outdir = f"{RESULTS}/homer",
        genome = lambda wildcards: wildcards.group.split("_")[-1]
    log:
        f"{LOGS}/homer/{{group}}_annotate.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: config["annotate_peaks"]["mem_mb"] * attempt,
        runtime = lambda wildcards, attempt: config["annotate_peaks"]["runtime"] * attempt,
    shell:
        '''
        exec > {log} 2>&1
        mkdir -p {params.outdir} 
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/configureHomer.pl -install {params.genome} 
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/annotatePeaks.pl {input} {params.genome} > {output}
        '''

rule homer_find_motifs_genome:
    input:
        f"{RESULTS}/homer/{{group}}_annotate.txt"
    output:
        multiext(f"{RESULTS}/homer/{{group}}/", "homerResults.html", "knownResults.html")
    conda:
        "../envs/data_analysis.yml"
    params:
        outdir = f"{RESULTS}/homer",
        genome = lambda wildcards: wildcards.group.split("_")[-1],
        args = config["find_motifs_genome"]["args"]
    threads:
        int(config["find_motifs_genome"]["threads"])
    log:
        f"{LOGS}/homer/{{group}}_findMotifs.log"
    resources:
        tmpdir=TEMP,
        cpus_per_task= lambda wildcards,threads: threads,
        mem_mb= lambda wildcards,attempt: config["find_motifs_genome"]["mem_mb"] * attempt,
        runtime= lambda wildcards,attempt: config["find_motifs_genome"]["runtime"] * attempt,
    shell:
        '''
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/findMotifsGenome.pl {input} {params.genome} {params.outdir}/{wildcards.group} -p {threads} {params.args}
        '''

rule plot_annotated_peaks:
    input:
        f"{RESULTS}/homer/{{sample}}_annotate.txt"
    output:
        distribution_plot = f"{RESULTS}/plots/{{sample}}_distribution.png",
        gene_distribution_plot =  f"{RESULTS}/plots/{{sample}}_genes.png"
    params:
        threshold_fraction = 0.01
    conda:
        "../envs/data_analysis.yml"
    log:
        f"{LOGS}/plot_annotated_peaks/{{sample}}.log"
    script:
        "../scripts/plot_annotated_peaks.py"
