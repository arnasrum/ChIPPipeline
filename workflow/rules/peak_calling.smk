import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import get_macs_input

for key, value in get_macs_input().items():
    rule:
        name: f"macs_callpeak_{key}"
        input:
            control = [*map(lambda file: f"results/{config['duplicate']}/" + file + ".bam", value["control"])],
            treatment = [*map(lambda file: f"results/{config['duplicate']}/" + file + ".bam",value["treatment"])],
        output:
            f"results/macs3/{key}_model.r",
            f"results/macs3/{key}_peaks.narrowPeak",
            f"results/macs3/{key}_peaks.xls",
            f"results/macs3/{key}_summits.bed"
        params:
            args = config["macs3"]["args"],
            name = key
        conda:
            "../envs/peak_calling.yml"
        shell:
            """
            macs3 callpeak -c {input.control} -t {input.treatment} --outdir results/macs3 --name {params.name} {params.args} 
            """

rule deeptools_bamCoverage:
    input:
        f"results/{config['duplicate']}/{{sample}}.bam"
    output:
        "results/deeptools/{sample}.bw"
    shell:
        """
        bamCovarage -b {input} -o {output}
        """