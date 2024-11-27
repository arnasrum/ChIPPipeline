import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import get_macs_input


macs_input = get_macs_input()
for sample, replicates in macs_input.items():
    for replicate, value in replicates.items():
        rule:
            name: f"peak_calling_macs3_callpeak_{sample}_rep{replicate}"
            input:
                control = [*map(lambda file: f"results/{config['duplicate_processor']}/" + file + ".bam", value["control"])],
                treatment = [*map(lambda file: f"results/{config['duplicate_processor']}/" + file + ".bam",value["treatment"])],
            output:
                multiext(f"results/macs3/{sample}_rep{replicate}","_peaks.xls","_summits.bed","_peaks.narrowPeak")
                if value["peak_type"] == "narrow" else
                multiext(f"results/macs3/{sample}_rep{replicate}_peaks",".xls",".broadPeak",".gappedPeak")

            params:
                args = config["macs3"]["args"],
                broad_peaks = value["peak_type"] == "broad",
                paired_end = config['paired_end'],
                name = f"{sample}_rep{replicate}"
            conda:
                "../envs/peak_calling.yml"
            shell:
                """
                inputOptions=''
                shopt -s nocasematch
                if [[ {params.broad_peaks} =~ true ]]; then
                    inputOptions+='--broad '
                fi 
                if [[ {params.paired_end} =~ true ]]; then
                    inputOptions+='-f BAMPE '
                else
                    inputOptions+='-f BAM '
                fi 
                macs3 callpeak -c {input.control} -t {input.treatment} --outdir results/macs3 --name {params.name} {params.args} $inputOptions
                """

rule deeptools_bamCoverage:
    input:
        bam = f"results/{config['duplicate_processor']}/{{sample}}.bam",
        bam_index = f"results/{config['duplicate_processor']}/{{sample}}.bam.bai"
    output:
        "results/deeptools-bamCoverage/{sample}.bw"
    conda:
        "../envs/peak_calling.yml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output}
        """

def get_consensus_peak_input(sample: str) -> list[str]:
    if len(macs_input[sample]) < 2: raise Exception("Cannot calculate consensus peak for sample with only one replicate.", sample)
    macs_extension = "_peaks.narrowPeak" if macs_input[sample][replicate]["peak_type"] == "narrow" else "_peaks.broadPeak"
    return [*map(lambda replicate: f"results/{config['peak_caller']}/{sample}_rep{replicate}{macs_extension}", macs_input[sample].keys())]

rule bedtools_consensus_peak:
    input:
        a = lambda wildcards: get_consensus_peak_input(wildcards.sample)[0],
        b = lambda wildcards: get_consensus_peak_input(wildcards.sample)[1:]
    output:
        "results/bedtools/{sample}.consensusPeak"
    conda:
        "../envs/peak_calling.yml"
    shell:
        '''
        bedtools intersect -a {input.a} -b {input.b} -wa > {output}
        '''


def get_compute_matrix_input(name: str, replicate: str) -> dict[str:list[str]]:
    macs_extension = "_peaks.narrowPeak" if macs_input[name][replicate]["peak_type"] == "narrow" else ["_peaks.broadPeak"]
    beds = [*map(lambda ext: f"results/{config['peak_caller']}/{name}_rep{replicate}{ext}", macs_extension)]
    bigwigs = [*map(lambda sample: f"results/deeptools-bamCoverage/{sample}.bw", macs_input[name][replicate]["treatment"])]
    return {"bigwigs": bigwigs, "beds": beds}

rule deeptools_computeMatrix:
    input:
        beds = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["beds"],
        bigwigs = lambda wildcards: get_compute_matrix_input(wildcards.sample, wildcards.replicate)["bigwigs"]
    output:
        "results/deeptools/{sample}_rep{replicate}_matrix.gz"
    conda:
        "../envs/peak_calling.yml"
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
        "../envs/peak_calling.yml"
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
        "../envs/peak_calling.yml"
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
        region = "chr1:10,000,000-11,000,000"
    conda:
        "../envs/peak_calling.yml"
    shell:
        ''' 
        make_tracks_file -f {input.beds} {input.bigwigs} -o {output.tracks}
        pyGenomeTracks --tracks {output.tracks} --region {params.region} --outFileName {output.plot}
        ''' 
rule annotate_peaks:
    input:
        peak ="results/bedtools/{sample}.consensusPeak",
    output:
        "results/homer/{sample}_annotate.txt"
    conda:
        "../envs/peak_calling.yml"
    params:
        genome = config['genome']
    shell:
        '''
        keyword='{params.genome}'
        grep_output=$(perl $CONDA_PREFIX/share/homer/configureHomer.pl -list | grep "+")
        if echo "$grep_output" | grep -q "$keyword"; then
            echo "Keyword '$keyword' found in the output."
        else
            perl $CONDA_PREFIX/share/homer/configureHomer.pl -install {params.genome} 
            echo "Keyword '$keyword' not found in the output."
        fi 
        perl $CONDA_PREFIX/share/homer/bin/annotatePeaks.pl {input.peak} {params.genome} > {output}
        '''

rule findMotifsGenome:
    input:
        "results/homer/{sample}_annotate.txt"
    output:
        multiext("results/homer/{sample}", "homerResults.html", "knownResults.html")
    conda:
        "../envs/peak_calling.yml"
    params:
        genome=config['genome'],
        size = 200
    shell:
        '''
        perl $CONDA_PREFIX/share/homer/bin/findMotifsGenome.pl {input} {params.genome} results/homer -size {params.size}
        '''
