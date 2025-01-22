import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import SampleFileScripts
from set_config_options import set_module_options, set_output_paths
from input_scripts import get_all_input, get_macs_input, symlink_input


set_module_options(config)
set_output_paths(config)
sfs = SampleFileScripts(config)
file_info = sfs.make_sample_info()

RESULTS = config['results_path']
RESOURCES = config['resources_path']
LOGS = config['logs_path']
BENCHMARKS = config['benchmarks_path']
TEMP = config['temp_path']
macs_input = get_macs_input(config["json_path"])


read_extention = ["_1", "_2"] if sfs.is_paired_end() else [""]
def get_trim_input(sample: str) -> list[str]:
    if sfs.is_paired_end():
        input = [RESOURCES + f"/reads/{sample}_1.fastq", RESOURCES + f"/reads/{sample}_2.fastq"]
    else:
        input = [RESOURCES + f"/reads/{sample}.fastq"]
    return input

def alignment_input(sample: str) -> list[str]:
    file_info = SampleFileScripts.get_file_info(config["json_path"])
    if sfs.is_paired_end():
        if sample in [*map(lambda accession: file_info['public'][accession]['file_name'], file_info["public"].keys())]:
            reads = [f"{RESULTS}/{config['trimmer']}/{sample}_1.fastq", f"{RESULTS}/{config['trimmer']}/{sample}_2.fastq"]
        else:
            reads = [f"{RESULTS}/{config['trimmer']}/{file_info['provided'][sample]['file_name']}_1.fastq",
                     f"{RESULTS}/{config['trimmer']}/{file_info['provided'][sample]['file_name']}_2.fastq"]
    else:
        if sample in [*map(lambda accession: file_info['public'][accession]['file_name'],file_info["public"].keys())]:
            reads = [f"{RESULTS}/{config['trimmer']}/{sample}.fastq"]
        else:
            reads = [f"{RESULTS}/{config['trimmer']}/{file_info['provided'][sample]['file_name']}.fastq"]
    return reads

def get_consensus_peak_input(sample: str) -> list[str]:
    peak_types = [*map(lambda replicate: macs_input[sample][replicate]["peak_type"],macs_input[sample])]
    if peak_types.count(
        peak_types[0]) != len(peak_types): raise ValueError(f"Peak types for {sample} do not match")
    macs_extension = "_peaks.narrowPeak" if peak_types[0] == "narrow" else "_peaks.broadPeak"
    return [*map(lambda replicate: f"{RESULTS}/{config['peak_caller']}/{sample}_rep{replicate}{macs_extension}",
        macs_input[sample].keys())]