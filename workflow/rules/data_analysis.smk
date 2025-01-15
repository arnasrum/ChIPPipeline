import sys
sys.path.append("workflow/scripts")
from input_scripts import get_macs_input

RESULTS: str = config['results_path']
LOGS: str = config['logs_path']
BENCHMARKS: str = config['benchmarks_path']
TEMP: str = config['temp_path']

macs_input = get_macs_input(config["json_path"])
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
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        bamCoverage -b {input.bam} -o {output}
        """

def get_consensus_peak_input(sample: str) -> list[str]:
    peak_types = [*map(lambda replicate: macs_input[sample][replicate]["peak_type"],macs_input[sample])]
    if peak_types.count(
        peak_types[0]) != len(peak_types): raise ValueError(f"Peak types for {sample} do not match")
    macs_extension = "_peaks.narrowPeak" if peak_types[0] == "narrow" else "_peaks.broadPeak"
    return [*map(lambda replicate: f"{RESULTS}/{config['peak_caller']}/{sample}_rep{replicate}{macs_extension}",
        macs_input[sample].keys())]

rule bedtools_consensus_peak:
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


def get_compute_matrix_input(name: str, replicate: str) -> dict[str:list[str]]:
    macs_extension = ["_peaks.narrowPeak"] if macs_input[name][replicate]["peak_type"] == "narrow" else ["_peaks.broadPeak"]
    beds = [*map(lambda ext: f"{RESULTS}/{config['peak_caller']}/{name}_rep{replicate}{ext}", macs_extension)]
    bigwigs = [*map(lambda sample: f"{RESULTS}/deeptools-bamCoverage/{sample}.bw", macs_input[name][replicate]["treatment"])]
    return {"bigwigs": bigwigs, "beds": beds}

rule deeptools_computeMatrix:
    input:
        beds = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["beds"],
        bigwigs = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["bigwigs"]
    output:
        RESULTS + "/deeptools/{sample}_rep{replicate}_matrix.gz"
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
        beds = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["beds"],
        bigwigs = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["bigwigs"]
    output:
        tracks = temp(RESULTS + "/pyGenomeTracks/{sample}_rep{replicate}_tracks.ini"),
        plot = RESULTS + "/pyGenomeTracks/{sample}_rep{replicate}.png"
    params:
        region = config["plot_regions"],
    conda:
        "../envs/data_analysis.yml"
    log:
        LOGS + "/pyGenomeTracks/{sample}_rep{replicate}.log"
    resources:
        tmpdir=TEMP
    shell:
        ''' 
        exec > {log} 2>&1
        make_tracks_file -f {input.beds} {input.bigwigs} -o {output.tracks}
        pyGenomeTracks --tracks {output.tracks} --region {params.region} --outFileName {output.plot}
        '''

rule annotate_peaks:
    input:
        lambda wildcards: f"{RESULTS}/bedtools/{wildcards.sample}.consensusPeak" if len(get_consensus_peak_input(wildcards.sample)) > 1 else get_consensus_peak_input(wildcards.sample)
    output:
        RESULTS + "/homer/{sample}_annotate.txt"
    conda:
        "../envs/data_analysis.yml"
    params:
        outdir = RESULTS + "/homer",
        genome = config['genome']
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

rule findMotifsGenome:
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
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/findMotifsGenome.pl ./{input} -p {threads} {params.genome} {params.outdir}/{wildcards.sample} -size {params.size}
        '''

