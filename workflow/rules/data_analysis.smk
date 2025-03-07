rule deeptools_bamCoverage:
    input:
        bam = RESULTS + "/" + config['duplicate_processor'] + "/{sample}.bam",
        bam_index = RESULTS + "/" + config['duplicate_processor'] + "/{sample}.bai",
    output:
        RESULTS + "/deeptools-bamCoverage/{sample}.bw"
    conda:
        "../envs/data_analysis.yml"
    log:
        LOGS + "/bamCoverage/{sample}.log"
    resources:
        tmpdir=TEMP
    threads:
        12
    shell:
        """
        exec > {log} 2>&1
        bamCoverage -p {threads} -b {input.bam} -o {output}
        """

rule deeptools_computeMatrix:
    input:
        beds = lambda wildcards:[f"{RESULTS}/{config['peak_caller']}/{wildcards.sample}.bed"],
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
        8
    resources:
        tmpdir=TEMP
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
        beds = lambda wildcards:[f"{RESULTS}/{config['peak_caller']}/{wildcards.sample}.bed"],
        bigwigs = lambda wildcards: f"{RESULTS}/deeptools-bamCoverage/{wildcards.sample}.bw"
    output:
        tracks = temp(RESULTS + "/pyGenomeTracks/{sample}_tracks.ini"),
        plot = RESULTS + "/pyGenomeTracks/{sample}.png"
    params:
        region = config["plot_regions"],
        args = config["pyGenomeTracks"]["args"]
    conda:
        "../envs/data_analysis.yml"
    log:
        LOGS + "/pyGenomeTracks/{sample}.log"
    resources:
        tmpdir=TEMP
    shell:
        ''' 
        exec > {log} 2>&1
        make_tracks_file -f {input.beds} {input.bigwigs} -o {output.tracks}
        python3 workflow/scripts/rename_tracks.py {output.tracks}
        pyGenomeTracks --tracks {output.tracks} --region {params.region} --outFileName {output.plot} {params.args}
        '''

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
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        bedtools intersect -a {input.a} -b {input.b} -wa > {output}
        '''

rule homer_annotate_peaks:
    input:
        lambda wildcards: f"{RESULTS}/bedtools/{wildcards.sample}.consensusPeak"
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
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        mkdir -p {params.outdir} 
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/configureHomer.pl -install {params.genome} 
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/annotatePeaks.pl {input} {params.genome} > {output}
        '''

rule homer_findMotifsGenome:
    input:
        RESULTS + "/homer/{sample}_annotate.txt"
    output:
        multiext(RESULTS + "/homer/{sample}/", "homerResults.html", "knownResults.html")
    conda:
        "../envs/data_analysis.yml"
    params:
        outdir = RESULTS + "/homer",
        genome = config['genome'],
        size = 200
    threads:
        12
    log:
        LOGS + "/homer/{sample}_findMotifs.log"
    resources:
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/findMotifsGenome.pl {input} {params.genome} {params.outdir}/{wildcards.sample} -size {params.size} -p {threads} 
        '''

rule plot_annotated_peaks:
    input:
        RESULTS + "/homer/{sample}_annotate.txt"
    output:
        distribution_plot = RESULTS + "/plots/{sample}_distribution.png",
        gene_distribution_plot =  RESULTS + "/plots/{sample}_genes.png"
    params:
        threshold_fraction = 0.01
    conda:
        "../envs/data_analysis.yml"
    script:
        "../scripts/plot_annotated_peaks.py"
