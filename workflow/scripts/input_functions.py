from workflow.scripts.pipeline_configuration import PipelineConfiguration

def trimmer_input(sample: str, resource_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    fastq_file_extensions = ["_1.fastq.gz", "_2.fastq.gz"] if pipeline_config.is_paired_end(sample) else [".fastq.gz"]
    return [f"{resource_path.rstrip('/')}/reads/{sample}{extension}" for extension in fastq_file_extensions]

def alignment_input(file_name: str, result_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    fastq_file_extensions = ["_1.fastq.gz", "_2.fastq.gz"] if pipeline_config.is_paired_end(file_name) else [".fastq.gz"]
    return [f"{result_path}/{pipeline_config.get_config_option('trimmer')}/{file_name}{extension}" for extension in fastq_file_extensions]

def peak_calling_input(sample: str, result_path: str, pipeline_config: PipelineConfiguration) -> dict[str, list[str]]:
    treatment_groups = pipeline_config.group_treatment_files()
    groups = [treatment_groups[pool][replicate]
        for pool in treatment_groups
        for replicate in treatment_groups[pool]
    ]
    treatment_group = next(filter(lambda group: sample in group, groups), None)
    path = f"{result_path}/{pipeline_config.get_config_option('duplicate_processor')}"
    return {"control": [f"{path}/{file}.bam" for file in pipeline_config.get_control_files(sample)],
            "treatment": [f"{path}/{file}.bam" for file in treatment_group]}

def get_consensus_peak_input(treatment_group: str, result_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    treatment_files = []
    treatment_groups = pipeline_config.group_treatment_files()
    for files in treatment_groups[treatment_group].values():
        treatment_files += [f"{result_path}/{pipeline_config.get_config_option('peak_caller')}/{files[0]}_peaks.narrowPeak"]

    if len(treatment_files) == 1:
        treatment_files.append(treatment_files[0])
    return treatment_files

def symlink_input(file_name: str, pipeline_config: PipelineConfiguration) -> str | None:
    pools = pipeline_config.sample_info["provided"].values()
    for pool in pools:
        for sample in pool:
            if sample['file_name'] == file_name:
                return sample
    return None

def concatenate_runs_input(file_name: str, resource_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    extensions = ["_1", "_2"] if pipeline_config.is_paired_end(file_name) else [""]
    runs = pipeline_config.get_sample_entry_by_file_name(file_name)["runs"]
    return [
        f"{resource_path}/reads/{run}{extension}.fastq"
        for extension in extensions
        for run in runs
    ]

def join_read_files(runs: list, paired_end: bool, resource_path: str) -> tuple[str] | list[str]:
    if paired_end:
        return [" ".join([f"{resource_path}/reads/{run}_1.fastq" for run in runs]),
                " ".join([f"{resource_path}/reads/{run}_2.fastq" for run in runs])
        ]
    return " ".join([f"{resource_path}/reads/{run}.fastq" for run in runs]),

def get_all_input(result_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    input_files = []
    path = result_path
    treatment_groups = pipeline_config.group_treatment_files()
    all_files = pipeline_config.get_all_file_names()
    if pipeline_config.get_config_option_bool('generate_fastqc_reports'):
        input_files += [get_fastqc_output(file_name, result_path, pipeline_config) for file_name in all_files]
    for treatment_file in pipeline_config.get_treatment_files():
        if pipeline_config.get_sample_entry_by_file_name(treatment_file)["peak_type"] == "narrow":
            input_files += [f"{path}/pyGenomeTracks/{treatment_file}.png"]
        if pipeline_config.get_sample_entry_by_file_name(treatment_file)["peak_type"] == "broad":
            input_files += [f"{path}/{pipeline_config.get_config_option('peak_caller')}/{treatment_file}_peaks.broadPeak"]
        input_files.append(f"{path}/deeptools/{treatment_file}_heatmap.png")
        input_files.append(f"{path}/deeptools/{treatment_file}_profile.png")
    for group, replicates in treatment_groups.items():
        allow_append = all(
            pipeline_config.get_sample_entry_by_file_name(sample)["peak_type"] == "narrow"
            for replicate in replicates
            for sample in replicates[replicate]
        )
        if allow_append:
            input_files.append(f"{path}/homer/{group}/homerResults.html")
            input_files.append(f"{path}/plots/{group}_genes.png")
    return input_files

def get_fastqc_output(file_name: str, result_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    base_paths = [f"{result_path}/fastqc/unprocessed/{file_name}",
                  f"{result_path}/fastqc/{pipeline_config.get_config_option('trimmer')}/{file_name}"
    ]
    single_end_extensions = ["_fastqc.zip", "_fastqc.html"]
    paired_end_extensions = ["_1_fastqc.zip", "_2_fastqc.zip", "_1_fastqc.html", "_2_fastqc.html"]
    extensions = paired_end_extensions if pipeline_config.is_paired_end(file_name) else single_end_extensions
    output_files = [
        f"{path}{extension}"
        for path in base_paths
        for extension in extensions
    ]
    return output_files

def flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict
