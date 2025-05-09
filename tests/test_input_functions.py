import pytest
from unittest.mock import MagicMock
from workflow.scripts.input_functions import *

@pytest.fixture
def mock_pipeline_config():
    config = MagicMock()
    # Define any necessary mock configurations or behaviors here
    return config


def test_trimmer_input_base_case_paired_end(mock_pipeline_config):
    mock_pipeline_config.is_paired_end.return_value = True
    sample = "sample1"
    resource_path = "resources"

    expected_outputs = [
        "resources/reads/sample1_1.fastq.gz",
        "resources/reads/sample1_2.fastq.gz"
    ]
    assert trimmer_input(sample, resource_path, mock_pipeline_config) == expected_outputs

    sample = "sample2"
    expected_outputs = [
        "resources/reads/sample2_1.fastq.gz",
        "resources/reads/sample2_2.fastq.gz"
    ]
    assert trimmer_input(sample, resource_path, mock_pipeline_config) == expected_outputs

def test_trimmer_input_new_path_paired_end(mock_pipeline_config):
    mock_pipeline_config.is_paired_end.return_value = True
    file_name = "sample"
    resource_path = "new_path/"
    expected_outputs = [
        "new_path/reads/sample_1.fastq.gz",
        "new_path/reads/sample_2.fastq.gz"
    ]
    assert trimmer_input(file_name, resource_path, mock_pipeline_config) == expected_outputs

def test_trimmer_input_base_case_single_end(mock_pipeline_config):
    mock_pipeline_config.is_paired_end.return_value = False
    file_name = "sample"
    resource_path = "results/"
    expected_outputs = [
        "results/reads/sample.fastq.gz",
    ]
    assert trimmer_input(file_name, resource_path, mock_pipeline_config) == expected_outputs

def test_alignment_input_base_case_paired_end(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "cutadapt"
    mock_pipeline_config.is_paired_end.return_value = True
    file_name = "file_name_A"
    results_path = "results"
    expected_outputs = [
        "results/cutadapt/file_name_A_1.fastq.gz",
        "results/cutadapt/file_name_A_2.fastq.gz"
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

def test_alignment_input_new_file_name_paired_end(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "cutadapt"
    mock_pipeline_config.is_paired_end.return_value = True
    results_path = "results"
    file_name = "file_name_A"

    expected_outputs = [
        "results/cutadapt/file_name_A_1.fastq.gz",
        "results/cutadapt/file_name_A_2.fastq.gz"
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

    file_name = "file_name_B"
    expected_outputs = [
        "results/cutadapt/file_name_B_1.fastq.gz",
        "results/cutadapt/file_name_B_2.fastq.gz"
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

def test_alignment_input_new_path_paired_end(mock_pipeline_config):
    mock_pipeline_config.is_paired_end.return_value = True
    mock_pipeline_config.get_config_option.return_value = "cutadapt"
    file_name = "file_name"
    results_path = "results"
    expected_outputs = [
        "results/cutadapt/file_name_1.fastq.gz",
        "results/cutadapt/file_name_2.fastq.gz"
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

    results_path = "new_path"
    expected_outputs = [
        "new_path/cutadapt/file_name_1.fastq.gz",
        "new_path/cutadapt/file_name_2.fastq.gz"
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

def test_alignment_input_new_trimmer_paired_end(mock_pipeline_config):
    mock_pipeline_config.is_paired_end.return_value = True
    mock_pipeline_config.get_config_option.return_value = "fastp"
    file_name = "file_name"
    results_path = "results"
    expected_outputs = [
        "results/fastp/file_name_1.fastq.gz",
        "results/fastp/file_name_2.fastq.gz"
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

    mock_pipeline_config.get_config_option.return_value = "trimmomatic"
    expected_outputs = [
        "results/trimmomatic/file_name_1.fastq.gz",
        "results/trimmomatic/file_name_2.fastq.gz"
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

def test_alignment_input_base_case_single_end(mock_pipeline_config):
    mock_pipeline_config.is_paired_end.return_value = False
    mock_pipeline_config.get_config_option.return_value = "cutadapt"
    file_name = "file_name"
    results_path = "results"
    expected_outputs = [
        "results/cutadapt/file_name.fastq.gz",
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

def test_peak_calling_input_base_case(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "mark_duplicates"
    mock_pipeline_config.group_treatment_files.return_value = {
        'antibody_sample_genome' : {
            1: [
                'GSM0000000_antibody_sample_treatment_genome1_rep1'
            ],
            2: [
                'GSM0000000_antibody_sample_treatment_genome2_rep2'
            ],
        }
    }
    mock_pipeline_config.get_control_files.side_effect = lambda key: {
        'GSM0000000_antibody_sample_treatment_genome1_rep1': ['control1'],
        'GSM0000000_antibody_sample_treatment_genome2_rep2': ['control2']
    }.get(key)

    result_path = "results"

    expected_outputs = {
            'control': [f"results/mark_duplicates/control1.bam"],
            'treatment': [f"results/mark_duplicates/GSM0000000_antibody_sample_treatment_genome1_rep1.bam"]
        }
    assert expected_outputs == peak_calling_input('GSM0000000_antibody_sample_treatment_genome1_rep1', result_path, mock_pipeline_config)
    
def test_peak_calling_input_treatment_pool(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "mark_duplicates"
    mock_pipeline_config.group_treatment_files.return_value = {
        'antibody_sample_genome' : {
            1: [
                'GSM0000000_antibody_sample_treatment_genome1_rep1',
                'GSM0000002_antibody_sample_treatment_genome1_rep1'
            ],
            2: [
                'GSM0000000_antibody_sample_treatment_genome2_rep2'
            ],
        }
    }
    mock_pipeline_config.get_control_files.side_effect = lambda key: {
        'GSM0000000_antibody_sample_treatment_genome1_rep1': ['control1'],
        'GSM0000000_antibody_sample_treatment_genome2_rep2': ['control2']
    }.get(key)

    result_path = "results"

    expected_outputs = {
            'control': ["results/mark_duplicates/control1.bam"],
            'treatment': ["results/mark_duplicates/GSM0000000_antibody_sample_treatment_genome1_rep1.bam",
                          f"results/mark_duplicates/GSM0000002_antibody_sample_treatment_genome1_rep1.bam"
            ]
        }
    assert expected_outputs == peak_calling_input('GSM0000000_antibody_sample_treatment_genome1_rep1', result_path, mock_pipeline_config)

def test_peak_calling_input_control_pool(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "mark_duplicates"
    mock_pipeline_config.group_treatment_files.return_value = {
        'antibody_sample_genome': {
            1: [
                'GSM0000000_antibody_sample_treatment_genome1_rep1',
            ],
            2: [
                'GSM0000000_antibody_sample_treatment_genome2_rep2'
            ],
        }
    }
    mock_pipeline_config.get_control_files.side_effect = lambda key: {
        'GSM0000000_antibody_sample_treatment_genome1_rep1': ['control1', 'control3'],
        'GSM0000000_antibody_sample_treatment_genome2_rep2': ['control2']
    }.get(key)
    result_path = "results"

    expected_outputs = {
        'control': ["results/mark_duplicates/control1.bam", "results/mark_duplicates/control3.bam"],
        'treatment': ["results/mark_duplicates/GSM0000000_antibody_sample_treatment_genome1_rep1.bam"]
    }
    assert expected_outputs == peak_calling_input('GSM0000000_antibody_sample_treatment_genome1_rep1', result_path, mock_pipeline_config)