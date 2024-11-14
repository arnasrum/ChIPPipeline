import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import get_macs_input

for key, value in get_macs_input().items():
    rule:
        name: f"macs_callpeak_{key[0]}_rep{key[1]}"
        input:
            control = [*map(lambda file: "results/picard-MarkDuplicates/" + file + ".bam", value["control"])],
            treatment = [*map(lambda file: "results/picard-MarkDuplicates/" + file + ".bam",value["treatment"])],
        output:
            f"results/macs3/{key[0]}_rep{key[1]}_model.r",
            f"results/macs3/{key[0]}_rep{key[1]}_peaks.narrowPeak",
            f"results/macs3/{key[0]}_rep{key[1]}_peaks.xls",
            f"results/macs3/{key[0]}_rep{key[1]}_summits.bed",
            "test.test"
        params:
            args = config["macs3"]["args"],
            name = f"{key[0]}_rep{key[1]}"
        conda:
            "../envs/peak_calling.yml"
        shell:
            """
            macs3 callpeak -c {input.control} -t {input.treatment} --outdir results/macs3 --name {params.name} {params.args} 
            """