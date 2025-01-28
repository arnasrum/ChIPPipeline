import sys
sys.path.append("workflow/scripts")
from pipeline_configuration import PipelineConfiguration
from set_config_options import set_module_options, set_output_paths
from input_scripts import get_macs_input, symlink_input, flatten_dict


set_module_options(config)
set_output_paths(config)
sfs = PipelineConfiguration(config)
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

def extract_files(sample, replicate, type) -> list[str]:
    return [*map(lambda file: f"{RESULTS}/{config['duplicate_processor']}/" + file + ".bam", get_macs_input(config['json_path'])[sample][replicate][type])]

def reference_genome_input():
    if os.path.isfile(config["genome"]):
        return config["genome"]
    return ""

def get_fastqc_output() -> list[str]:
    file_names = sfs.get_all_file_names()
    file_paths_unprocessed = [f"{RESULTS}/fastqc/unprocessed/{file_name}" for file_name in file_names]
    file_paths_trimmed = [f"{RESULTS}/fastqc/{config['trimmer']}/{file_name}" for file_name in file_names]

    single_end_extensions = ["_fastqc.zip", "_fastqc.html"]
    paired_end_extensions = ["_1_fastqc.zip", "_2_fastqc.zip", "_1_fastqc.html", "_2_fastqc.html"]
    extensions = paired_end_extensions if sfs.is_paired_end() else single_end_extensions

    output_files = [f"{path}{extension}"
                    for path in file_paths_unprocessed + file_paths_trimmed
                    for extension in extensions
    ]
    return output_files


def get_all_input(config):
    input_files = []
    path = config["results_path"]
    for key, value in macs_input.items():
        for replicate in value:
            input_files.append(f"{path}/deeptools/{key}_rep{replicate}_heatmap.png")
            input_files.append(f"{path}/deeptools/{key}_rep{replicate}_profile.png")
            input_files.append(f"{path}/pyGenomeTracks/{key}_rep{replicate}.png")

    # Gets narrow peak samples
    input_files += get_fastqc_output()
    input_files += [*map(lambda sample: f"{config['results_path']}/homer/{sample}/homerResults.html",
                         filter(lambda sample: macs_input[sample][next(iter(macs_input[sample].keys()))]["peak_type"] == "narrow", macs_input.keys()))]
    return input_files