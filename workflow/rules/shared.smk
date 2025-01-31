import sys
sys.path.append("workflow/scripts")
from pipeline_configuration import PipelineConfiguration
from set_config_options import set_module_options, set_output_paths

set_module_options(config)
set_output_paths(config)
sfs = PipelineConfiguration(config)
file_info = sfs.make_sample_info()

RESULTS = config['results_path']
RESOURCES = config['resources_path']
LOGS = config['logs_path']
BENCHMARKS = config['benchmarks_path']
TEMP = config['temp_path']
fastq_file_extensions = ["_1.fastq", "_2.fastq"] if sfs.is_paired_end() else [".fastq"]


def get_trim_input(sample: str) -> list[str]:
    return [f"{RESOURCES}/reads/{sample}{extension}" for extension in fastq_file_extensions]

def alignment_input(file_name: str) -> list[str]:
    return [f"{RESULTS}/{config["trimmer"]}/{file_name}{extension}" for extension in fastq_file_extensions]

def get_consensus_peak_input(sample: str) -> list[str]:
    peak_types = [*map(lambda replicate: macs_input[sample][replicate]["peak_type"],macs_input[sample])]
    if peak_types.count(
        peak_types[0]) != len(peak_types): raise ValueError(f"Peak types for {sample} do not match")
    macs_extension = "_peaks.narrowPeak" if peak_types[0] == "narrow" else "_peaks.broadPeak"
    return [*map(lambda replicate: f"{RESULTS}/{config['peak_caller']}/{sample}_rep{replicate}{macs_extension}",
        macs_input[sample].keys())]

def macs_input_func(sample, replicate, type) -> list[str]:
    return [*map(lambda file: f"{RESULTS}/{config['duplicate_processor']}/" + file + ".bam", get_macs_input()[sample][replicate][type])]

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
    input_files += [*map(lambda sample: f"{RESULTS}/homer/{sample}/homerResults.html",
                         filter(lambda sample: macs_input[sample][next(iter(macs_input[sample].keys()))]["peak_type"] == "narrow", macs_input.keys()))]
    return input_files


def symlink_input(json_path: str, file_name: str) -> None:
    with open(json_path) as file:
        samples = sfs.sample_info["provided"]
    return next((item for item in samples.values() if item["file_name"] == file_name), None)

def flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict

def get_macs_input() -> dict[str: dict]:
    """
    Groups input samples into replicates and treatment/control for sample and mark.
    Two unique samples with same type, sample, replicate and mark, get pooled together.
    """
    samples = sfs.sample_info
    macs_input = {}
    control_files = []

    # Process all sample entries
    for key, entry in flatten_dict(samples).items():
        replicate = str(entry["replicate"])
        sample_type = entry["type"]
        sample_mark = entry["mark"]
        sample_file = entry["file_name"]
        sample_name = entry["sample"]

        if sample_type == "treatment":
            file_key = f"{sample_mark}_{sample_name}"
            macs_input.setdefault(file_key,{}).setdefault(replicate,{
                "treatment": [],
                "control": [],
                "peak_type": entry['peak_type']
            })
            macs_input[file_key][replicate]["treatment"].append(sample_file)

        elif sample_type == "control":
            control_files.append((sample_name, replicate, sample_file))

        else:
            raise Exception(f"Entry type unrecognized for {sample_file}")

    # Associate control files with their respective treatments
    for sample_name, replicate, file_name in control_files:
        for file_key in filter(lambda k: sample_name in k,macs_input.keys()):
            macs_input[file_key][replicate]["control"].append(file_name)

    # Remove any replicates with missing treatment or control
    for file_key, replicates in macs_input.items():
        invalid_replicates = [
            replicate for replicate, data in replicates.items()
            if not data["treatment"] or not data["control"]
        ]
        for replicate in invalid_replicates:
            del replicates[replicate]

    return macs_input

macs_input = get_macs_input()
