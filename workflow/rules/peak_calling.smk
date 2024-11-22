import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import get_macs_input


macs_input = get_macs_input()
for key, value in macs_input.items():
    rule:
        name: f"macs_callpeak_{key}"
        input:
            control = [*map(lambda file: f"results/{config['duplicate_processor']}/" + file + ".bam", value["control"])],
            treatment = [*map(lambda file: f"results/{config['duplicate_processor']}/" + file + ".bam",value["treatment"])],
        output:
            multiext(f"results/macs3/{key}", "_peaks.xls", "_model.r", "_summits.bed", "_peaks.narrowPeak")
            if value["peak_type"] == "narrow" else
            multiext(f"results/macs3/{key}_peaks", ".xls", ".broadPeak", ".gappedPeak")

        params:
            args = config["macs3"]["args"],
            broad_peaks = value["peak_type"] == "broad",
            paired_end = config['paired_end'],
            name = key
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
        "results/deeptools/{sample}.bw"
    conda:
        "../envs/peak_calling.yml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output}
        """


def get_compute_matrix_input(name: str) -> dict[str:list[str]]:
    macs_extensions = ["_summits.bed","_peaks.narrowPeak"] if macs_input[name]["peak_type"] == "narrow" else ["_peaks.broadPeak","_peaks.gappedPeak"]
    bigwigs = [*map(lambda ext: f"results/{config['peak_caller']}/{name}{ext}", macs_extensions)]
    beds = [*map(lambda sample: f"results/{config['duplicate_processor']}/{sample}.bam", macs_input[name]["treatment"])]
    return {"bigwigs": bigwigs, "beds": beds}

rule deeptools_computeMatrix:
    input:
        beds = lambda wildcards: [bed for bed in get_compute_matrix_input(wildcards.sample)["beds"]],
        bigwigs = lambda wildcards: [bigwig for bigwig in get_compute_matrix_input(wildcards.sample)["bigwigs"]]
    output:
        "results/deeptools/{sample}_matrix.gz"
    conda:
        "../envs/peak_calling.yml"
    params:
        mode = "reference-point",
        args = ""
    shell:
        """
        computeMatrix {params.mode} -S {input.bigwigs} -R {input.beds} -o {output} {params.args}
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

