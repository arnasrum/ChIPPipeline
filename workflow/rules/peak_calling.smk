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
            multiext(f"results/macs3/{key}", "_peaks.xls", "_model.r", "_summits.bed", "_peaks.narrowPeak") if value["peak_type"]
            else multiext(f"results/macs3/{key}", "_peaks.xls", "_model.r", "_gapped.bed", "_peaks.broadPeak")

        params:
            args = config["macs3"]["args"],
            paired_end = config['paired_end'],
            name = key
        conda:
            "../envs/peak_calling.yml"
        shell:
            """
            if [[ {params.paired_end} ]]; then
                {params.args} += '-f BAMPE'
            else
                {params.args} += '-f BAM'
            fi
            macs3 callpeak -c {input.control} -t {input.treatment} --outdir results/macs3 --name {params.name} {params.args}
            """

rule deeptools_bamCoverage:
    input:
        bam = f"results/{config['aligner']}/{{sample}}.bam",
        bam_index = f"results/{config['aligner']}/{{sample}}.bam.bai"
    output:
        "results/deeptools/{sample}.bw"
    conda:
        "../envs/peak_calling.yml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output}
        """