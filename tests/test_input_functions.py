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



