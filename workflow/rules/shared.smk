from scripts.pipeline_configuration import PipelineConfiguration
from scripts.set_config_options import set_module_options, set_output_paths

set_module_options(config)
set_output_paths(config)
sfs = PipelineConfiguration(config)
file_info = sfs.make_sample_info()
genome = sfs.get_genome_code()
treatment_groups = sfs.group_treatment_files()

RESULTS = config['results_path']
RESOURCES = config['resources_path']
LOGS = config['logs_path']
BENCHMARKS = config['benchmarks_path']
TEMP = config['temp_path']
fastq_file_extensions = ["_1.fastq", "_2.fastq"] if sfs.is_paired_end() else [".fastq"]

def trimmer_input(sample: str) -> list[str]:
    return [f"{RESOURCES}/reads/{sample}{extension}" for extension in fastq_file_extensions]

def alignment_input(file_name: str) -> list[str]:
    return [f"{RESULTS}/{config["trimmer"]}/{file_name}{extension}" for extension in fastq_file_extensions]

def get_consensus_peak_input(treatment_group: str) -> list[str]:
    treatment_files = []
    for files in treatment_groups[treatment_group].values():
        treatment_files += [f"{RESULTS}/{config['peak_caller']}/{file}.bed" for file in files]
    if len(treatment_files) == 1:
        treatment_files.append(treatment_files[0])
    return treatment_files

def macs_input(sample: str) -> dict[str, list[str]]:
    groups = [treatment_groups[pool][replicate]
        for pool in treatment_groups
        for replicate in treatment_groups[pool]
    ]
    treatment_group = next(filter(lambda group: sample in group, groups))
    path = f"{RESULTS}/{config['duplicate_processor']}"
    #print("control", [f"{path}/{file}.bam" for file in sfs.get_control_files(sample)])
    #print("treatment", [f"{path}/{file}.bam" for file in treatment_group])
    return {"control": [f"{path}/{file}.bam" for file in sfs.get_control_files(sample)],
            "treatment": [f"{path}/{file}.bam" for file in treatment_group]}

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
    input_files += get_fastqc_output()
    for treatment_file in sfs.get_treatment_files():
        input_files += [f"{path}/pyGenomeTracks/{treatment_file}.png"]
        input_files.append(f"{path}/deeptools/{treatment_file}_heatmap.png")
        input_files.append(f"{path}/deeptools/{treatment_file}_profile.png")
    for group, replicates in treatment_groups.items():
        allow_append = all(
            sfs.get_sample_entry_by_file_name(sample)["peak_type"] == "narrow"
            for replicate in replicates
            for sample in replicates[replicate]
        )
        if allow_append:
            input_files.append(f"{RESULTS}/homer/{group}/homerResults.html")
            input_files.append(f"{path}/plots/{group}_genes.png")
    return input_files

def get_macs_input(config):
    return [f"{RESULTS}/macs3/{treatment_file}.bed" for treatment_file in sfs.get_treatment_files()]



def symlink_input(json_path: str, file_name: str) -> None:
    with open(json_path) as file:
        samples = sfs.sample_info["provided"]
    return next((item for item in samples.values() if item["file_name"] == file_name), None)

def flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict
