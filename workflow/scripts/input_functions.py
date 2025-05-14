from workflow.scripts.pipeline_configuration import PipelineConfiguration
import os

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
                if pipeline_config.is_paired_end(file_name):
                    return[sample['read1']['path'], sample['read2']['path']]
                else:
                    return[sample['read1']['path']]
    return None

def concatenate_runs_input(file_name: str, resource_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    extensions = ["_1", "_2"] if pipeline_config.is_paired_end(file_name) else [""]
    runs = pipeline_config.get_sample_entry_by_file_name(file_name)["runs"]
    return [
        f"{resource_path}/reads/{run}{extension}.fastq"
        for extension in extensions
        for run in runs
    ]

def join_read_files(runs: list, paired_end: bool, resource_path: str) -> str | list[str]:
    resource_path = resource_path.rstrip('/')
    if paired_end:
        return [" ".join([f"{resource_path}/reads/{run}_1.fastq" for run in runs]),
                " ".join([f"{resource_path}/reads/{run}_2.fastq" for run in runs])
        ]
    return " ".join([f"{resource_path}/reads/{run}.fastq" for run in runs])


HOMER_SUPPORTED_GENOMES = ['panTro4', 'gorGor5', 'rn4', 'corn.AGP', 'ce6', 'mm10', 'hg19', 'galGal4', 'danRer10', 'ce10', 'panPan3', 'ci3', 'xenLae2', 'rn7', 'rheMac2', 'ce11', 'dm3', 'patens.ASM242', 'rice.IRGSP-1.0', 'hg17', 'taeGut2', 'anoGam1', 'sacCer3', 'petMar2', 'xenTro7', 'susScr11', 'rn6', 'apiMel2', 'AGP', 'hg18', 'panTro6', 'xenTro3', 'fr3', 'canFam5', 'tair10', 'rheMac10', 'galGal6', 'panTro3', 'dm6', 'aplCal1', 'papAnu2', 'susScr3', 'mm9', 'tetNig2', 'panPan2', 'papAnu4', 'hg38', 'gorGor6', 'xenTro10', 'rheMac3', 'gorGor3', 'mm39', 'anoGam3', 'mm8', 'gorGor4', 'panTro5', 'canFam6', 'sacCer2', 'canFam3', 'petMar3', 'galGal5', 'rn5', 'ci2', 'apiMel3', 'xenTro9', 'panPan1', 'danRer11', 'xenTro2', 'danRer7', 'strPur2', 'rheMac8']
def get_all_input(result_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    input_files = []
    path = result_path
    treatment_groups = pipeline_config.group_treatment_files()
    all_files = pipeline_config.get_all_file_names()
    if pipeline_config.get_config_option_bool('generate_fastqc_reports'):
        input_files += [get_fastqc_output(file_name, result_path, pipeline_config) for file_name in all_files]
    for treatment_file in pipeline_config.get_treatment_files():
        if pipeline_config.get_sample_entry_by_file_name(treatment_file)["peak_type"] == "narrow":
            input_files.append(f"{path}/pyGenomeTracks/{treatment_file}.png")
            input_files.append(f"{path}/deeptools/{treatment_file}_heatmap.png")
            input_files.append(f"{path}/deeptools/{treatment_file}_profile.png")
        if pipeline_config.get_sample_entry_by_file_name(treatment_file)["peak_type"] == "broad":
            input_files.append(f"{path}/{pipeline_config.get_config_option('peak_caller')}/{treatment_file}_peaks.broadPeak")

    for group, replicates in treatment_groups.items():
        allow_append = all(
            pipeline_config.get_sample_entry_by_file_name(sample)["peak_type"] == "narrow"
                and pipeline_config.get_sample_genome(sample) in HOMER_SUPPORTED_GENOMES
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

def get_pooled_treatment_samples(file_name: str, pipeline_config: PipelineConfiguration) -> list[str]:
    groups = pipeline_config.group_treatment_files()
    for group in groups:
        for replicate in groups[group]:
            if file_name in groups[group][replicate]:
                return groups[group][replicate]
    raise ValueError(f"No samples found for {file_name}")


def get_genome_path(genome: str, pipeline_config: PipelineConfiguration) -> str | None:
    genomes = [
        sample['genome']
        for pool in flatten_dict(pipeline_config.sample_info).values()
        for sample in pool
    ]
    for genome in genomes:
        if os.path.isfile(genome):
            return genome
    return None

def bigwig_compare_input(sample: str, results_path: str, pipeline_config: PipelineConfiguration) -> list[str]:
    base_path = f"{results_path}/deeptools-bamCoverage"
    input_files = {"treatment": [], "control": []}
    control_files = pipeline_config.get_control_files(sample)
    input_files['treatment'] = [f"{base_path}/{file_name}.bw" for file_name in get_pooled_treatment_samples(sample, pipeline_config)]
    if control_files:
        input_files["control"] = [f"{base_path}/{file_name}.bw" for file_name in control_files]
    return input_files

def flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict
