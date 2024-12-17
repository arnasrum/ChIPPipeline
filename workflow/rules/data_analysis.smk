import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import get_macs_input

macs_input = get_macs_input()
rule deeptools_bamCoverage:
    input:
        bam = f"results/{config['duplicate_processor']}/{{sample}}.bam",
        bam_index = f"results/{config['duplicate_processor']}/{{sample}}.bam.bai"
    output:
        "results/deeptools-bamCoverage/{sample}.bw"
    conda:
        "../envs/data_analysis.yml"
    log:
        "logs/bamCoverage/{sample}.log"
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
    return [*map(lambda replicate: f"results/{config['peak_caller']}/{sample}_rep{replicate}{macs_extension}",
        macs_input[sample].keys())]

rule bedtools_consensus_peak:
    input:
        a = lambda wildcards: get_consensus_peak_input(wildcards.sample)[0],
        b = lambda wildcards: get_consensus_peak_input(wildcards.sample)[1:]
    output:
        "results/bedtools/{sample}.consensusPeak"
    conda:
        "../envs/peak_calling.yml"
    log:
        "logs/bedtools-intersect/{sample}.log"
    shell:
        '''
        exec > {log} 2>&1
        bedtools intersect -a {input.a} -b {input.b} -wa > {output}
        '''


def get_compute_matrix_input(name: str, replicate: str) -> dict[str:list[str]]:
    macs_extension = ["_peaks.narrowPeak"] if macs_input[name][replicate]["peak_type"] == "narrow" else ["_peaks.broadPeak"]
    beds = [*map(lambda ext: f"results/{config['peak_caller']}/{name}_rep{replicate}{ext}", macs_extension)]
    bigwigs = [*map(lambda sample: f"results/deeptools-bamCoverage/{sample}.bw", macs_input[name][replicate]["treatment"])]
    return {"bigwigs": bigwigs, "beds": beds}

rule deeptools_computeMatrix:
    input:
        beds = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["beds"],
        bigwigs = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["bigwigs"]
    output:
        "results/deeptools/{sample}_rep{replicate}_matrix.gz"
    wildcard_constraints:
        replicate = r"[0-9]"
    conda:
        "../envs/data_analysis.yml"
    params:
        mode = "reference-point",
        args = config["computeMatrix"]["args"]
    threads:
        8
    shell:
        """
        mkdir -p results/deeptools
        computeMatrix {params.mode} -p {threads} -S {input.bigwigs} -R {input.beds} -o {output} {params.args}
        """

rule deeptools_plotHeatMap:
    input:
        "results/deeptools/{sample}_matrix.gz"
    output:
        "results/deeptools/{sample}_heatmap.png"
    conda:
        "../envs/data_analysis.yml"
    shell:
        """
        plotHeatmap -m {input} -o {output}
        """
rule deeptools_plotProfile:
    input:
        "results/deeptools/{sample}_matrix.gz"
    output:
        "results/deeptools/{sample}_profile.png"
    conda:
        "../envs/data_analysis.yml"
    shell:
        """
        plotProfile -m {input} -o {output}
        """

rule plot_genome_track:
    input:
        beds = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["beds"],
        bigwigs = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["bigwigs"]
    output:
        tracks = temp("results/pyGenomeTracks/{sample}_rep{replicate}_tracks.ini"),
        plot = "results/pyGenomeTracks/{sample}_rep{replicate}.png"
    params:
        region = "chr1:10,000,000-11,000,000",
    conda:
        "../envs/peak_calling.yml"
    log:
        "logs/pyGenomeTracks/{sample}_rep{replicate}.log"
    shell:
        ''' 
        exec > {log} 2>&1
        python3 workflow/scripts/rename_peaks.py
        make_tracks_file -f {input.beds} {input.bigwigs} -o {output.tracks}
        pyGenomeTracks --tracks {output.tracks} --region {params.region} --outFileName {output.plot}
        '''
rule annotate_peaks:
    input:
        peak ="results/bedtools/{sample}.consensusPeak",
    output:
        "results/homer/{sample}_annotate.txt"
    conda:
        "../envs/homer.yml"
    params:
        genome = config['genome']
    log:
        "logs/homer/{sample}_annotate.log"
    shell:
        '''
        exec > {log} 2>&1
        keyword='{params.genome}'
        grep_output=$(perl $CONDA_PREFIX/share/homer/configureHomer.pl -list | grep "+")
        if echo "$grep_output" | grep -q "$keyword"; then
            echo "Keyword '$keyword' found in the output."
        else
            perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/configureHomer.pl -install {params.genome} 
            echo "Keyword '$keyword' not found in the output."
        fi 
        mkdir -p results/homer
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/annotatePeaks.pl {input.peak} {params.genome} > {output}
        '''

rule findMotifsGenome:
    input:
        "results/homer/{sample}_annotate.txt"
    output:
        multiext("results/homer/{sample}/", "homerResults.html", "knownResults.html")
    conda:
        "../envs/homer.yml"
    params:
        genome=config['genome'],
        size = 200
    log:
        "logs/homer/{sample}_findMotifs.log"
    shell:
        '''
        exec > {log} 2>&1
        mkdir -p results/homer
        perl -I $CONDA_PREFIX/share/homer/bin $CONDA_PREFIX/share/homer/bin/findMotifsGenome.pl ./{input} {params.genome} results/homer/{wildcards.sample} -size {params.size}
        '''

