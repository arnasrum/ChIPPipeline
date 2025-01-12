import sys
sys.path.append("workflow/scripts")
from input_scripts import get_macs_input

RESULTS: str = config['results_path']
LOGS: str = config['logs_path']
TEMP: str = config['temp_path']
BENCHMARKS: str = config['benchmarks_path']

macs_input = get_macs_input(config["json_path"])
for sample, replicates in macs_input.items():
    for replicate, value in replicates.items():
        rule:
            name: f"peak_calling_macs3_callpeak_{sample}_rep{replicate}"
            input:
                control = [*map(lambda file: f"{RESULTS}/{config['duplicate_processor']}/" + file + ".bam", value["control"])],
                treatment = [*map(lambda file: f"{RESULTS}/{config['duplicate_processor']}/" + file + ".bam",value["treatment"])],
            output:
                multiext(f"{RESULTS}/macs3/{sample}_rep{replicate}","_peaks.xls","_summits.bed","_peaks.narrowPeak")
                if value["peak_type"] == "narrow" else
                multiext(f"{RESULTS}/macs3/{sample}_rep{replicate}_peaks",".xls",".broadPeak",".gappedPeak")

            params:
                args = config["macs3"]["args"],
                broad_peaks = value["peak_type"] == "broad",
                paired_end = config['paired_end'],
                peak_type = value["peak_type"],
                name = f"{sample}_rep{replicate}",
                outdir = f"{RESULTS}/macs3"

            conda:
                "../envs/peak_calling.yml"
            log:
                f"{LOGS}/macs3/{sample}_rep{replicate}.log"
            benchmark:
                f"{BENCHMARKS}/macs3/{sample}_rep{replicate}.log"
            resources:
                tmpdir=TEMP
            shell:
                """
                exec > {log} 2>&1
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
                macs3 callpeak --tempdir {resources.tmpdir} -c {input.control} -t {input.treatment} --outdir {params.outdir} --name {params.name} {params.args} $inputOptions
                python3 workflow/scripts/rename_peaks.py {params.outdir}/{params.name}_peaks.{params.peak_type}Peak
                """
