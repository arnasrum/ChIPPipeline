import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import get_macs_input

for key, value in get_macs_input().items():
    rule:
        name: f"macs_callpeak_{key}"
        input:
            control = [*map(lambda file: f"results/{config['duplicate_processor']}/" + file + ".bam", value["control"])],
            treatment = [*map(lambda file: f"results/{config['duplicate_processor']}/" + file + ".bam",value["treatment"])],
        output:
            multiext(f"results/macs3/{key}", "_peaks.xls", "_model.r", "_summits.bed", "_peaks.narrowPeak") if value["peak_type"] == "narrow"
            else multiext(f"results/macs3/{key}_peaks", ".xls", ".broadPeak", ".gappedPeak")

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

rule deeptools_computeMatrix:
    input:
        ""
    output:
        ""
    conda:
        "../envs/peak_calling.yml"
    shell:
        """
        computeMatrix 
        """


