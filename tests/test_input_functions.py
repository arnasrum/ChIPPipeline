import pytest
from unittest.mock import MagicMock
from workflow.scripts.input_functions import *

@pytest.fixture
def mock_pipeline_config():
    config = MagicMock()
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

def test_alignment_input_new_file_name_single_end(mock_pipeline_config):
    mock_pipeline_config.is_paired_end.return_value = False
    mock_pipeline_config.get_config_option.return_value = "cutadapt"
    file_name = "file_name"
    results_path = "results"
    expected_outputs = [
        "results/cutadapt/file_name.fastq.gz",
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs
    file_name = "new_file_name"
    expected_outputs = [
        "results/cutadapt/new_file_name.fastq.gz",
    ]
    assert alignment_input(file_name, results_path, mock_pipeline_config) == expected_outputs

def test_peak_calling_input_base_case(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "mark_duplicates"
    mock_pipeline_config.group_treatment_files.return_value = {
        "mark1_sample1_mouse-genome": {"1": ["file_name1"], "2": ["file_name2"]},
        "mark2_sample1_mouse-genome": {"1": ["file_name3"], "2": ["file_name4"]},
    }
    mock_pipeline_config.get_control_files.return_value = ["control_file1"]
    file_name = "file_name1"
    results_path = "results"
    expected_outputs = {
            "treatment": ["results/mark_duplicates/file_name1.bam"],
            "control": ["results/mark_duplicates/control_file1.bam"],
    }
    assert expected_outputs == peak_calling_input(file_name, results_path, mock_pipeline_config)

def test_peak_calling_input_treatment_pooled(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "mark_duplicates"
    mock_pipeline_config.group_treatment_files.return_value = {
        "mark1_sample1_mouse-genome": {"1": ["file_name1", "file_name5"], "2": ["file_name2"]},
        "mark2_sample1_mouse-genome": {"1": ["file_name3"], "2": ["file_name4"]},
    }
    mock_pipeline_config.get_control_files.return_value = ["control_file1"]
    file_name = "file_name1"
    results_path = "results"
    expected_outputs = {
            "treatment": ["results/mark_duplicates/file_name1.bam", "results/mark_duplicates/file_name5.bam"],
            "control": ["results/mark_duplicates/control_file1.bam"],
    }
    assert expected_outputs == peak_calling_input(file_name, results_path, mock_pipeline_config)

def test_peak_calling_input_control_pooled(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "mark_duplicates"
    mock_pipeline_config.group_treatment_files.return_value = {
        "mark1_sample1_mouse-genome": {"1": ["file_name1"], "2": ["file_name2"]},
        "mark2_sample1_mouse-genome": {"1": ["file_name3"], "2": ["file_name4"]},
    }
    mock_pipeline_config.get_control_files.return_value = ["control_file1", "control_file2"]
    file_name = "file_name1"
    results_path = "results"
    expected_outputs = {
            "treatment": ["results/mark_duplicates/file_name1.bam"],
            "control": ["results/mark_duplicates/control_file1.bam", "results/mark_duplicates/control_file2.bam"],
    }
    assert expected_outputs == peak_calling_input(file_name, results_path, mock_pipeline_config)

def test_consensus_peak_input_base_case(mock_pipeline_config):
    mock_pipeline_config.group_treatment_files.return_value = {
        "mark1_sample1_mouse-genome": {"1": ["file_name1"], "2": ["file_name2"]},
    }
    mock_pipeline_config.get_config_option.return_value = "peak_caller"
    results_path = "results"
    group = "mark1_sample1_mouse-genome"

    expected_outputs = [
        "results/peak_caller/file_name1_peaks.narrowPeak",
        "results/peak_caller/file_name2_peaks.narrowPeak"
    ]
    assert expected_outputs == get_consensus_peak_input(group, results_path, mock_pipeline_config)

def test_consensus_peak_input_single_replicate(mock_pipeline_config):
    mock_pipeline_config.group_treatment_files.return_value = {
        "mark1_sample1_mouse-genome": {"1": ["file_name1"]},
    }
    mock_pipeline_config.get_config_option.return_value = "peak_caller"
    results_path = "results"
    group = "mark1_sample1_mouse-genome"

    expected_outputs = [
        "results/peak_caller/file_name1_peaks.narrowPeak",
        "results/peak_caller/file_name1_peaks.narrowPeak"
    ]
    assert expected_outputs == get_consensus_peak_input(group, results_path, mock_pipeline_config)

def test_handle_provided_input_base_case(mock_pipeline_config):
    mock_pipeline_config.sample_info = {
        "provided": {
            "antibody_sample1_rep1_mm10_treatment": [
                {
                    "read1": {
                        "path": "data/treatment_sample_1.fq.gz",
                        "file_extension": ".fq.gz",
                        "file_name": "treatment_sample_1"
                    },
                    "read2": {
                        "path": "data/treatment_sample_2.fq.gz",
                        "file_extension": ".fq.gz",
                        "file_name": "treatment_sample_2"
                    },
                    "file_name": "antibody_sample1_treatment_rep1_mm10",
                    "type": "treatment",
                    "sample": "sample1",
                    "replicate": 1,
                    "mark": "antibody",
                    "peak_type": "narrow",
                    "genome": "data/mm10.fa.gz",
                    "paired_end": "true"
                }
            ],
            "sample1_rep1_mm10_control": [
                {
                    "read1": {
                        "path": "data/control_sample_1.fq.gz",
                        "file_extension": ".fq.gz",
                        "file_name": "control_sample_1"
                    },
                    "read2": {
                        "path": "data/control_sample_2.fq.gz",
                        "file_extension": ".fq.gz",
                        "file_name": "control_sample_2"
                    },
                    "file_name": "GSM0123450_sample1_control_rep1_mm10",
                    "type": "control",
                    "sample": "sample1",
                    "replicate": 1,
                    "mark": "",
                    "peak_type": "",
                    "genome": "mm10",
                    "paired_end": "true"
                }
            ],
        }
    }
    expected_output = ["data/treatment_sample_1.fq.gz", "data/treatment_sample_2.fq.gz"] 
    assert expected_output == handle_provided_input("antibody_sample1_treatment_rep1_mm10", mock_pipeline_config)

def test_handle_provided_input_no_results(mock_pipeline_config):
    mock_pipeline_config.sample_info = {"provided": {}}
    assert None ==  handle_provided_input("antibody_sample1", mock_pipeline_config)

def test_join_read_files_paired_end():
    reads = ["SRA111", "SRA222", "SRA333"]
    paired_end = True
    resource_path = "resources"

    expected_output = ["resources/reads/SRA111_1.fastq resources/reads/SRA222_1.fastq resources/reads/SRA333_1.fastq",
                       "resources/reads/SRA111_2.fastq resources/reads/SRA222_2.fastq resources/reads/SRA333_2.fastq"]

    assert expected_output == join_read_files(reads, paired_end, resource_path)

def test_join_read_files_single_end():
    reads = ["SRA111", "SRA222", "SRA333"]
    paired_end = False
    resource_path = "resources"

    expected_output = "resources/reads/SRA111.fastq resources/reads/SRA222.fastq resources/reads/SRA333.fastq"

    assert expected_output == join_read_files(reads, paired_end, resource_path)

def test_get_fastqc_output_paired_end(mock_pipeline_config):
        mock_pipeline_config.get_config_option.return_value = "trimmer"
        mock_pipeline_config.is_paired_end.return_value = True
        file_name = "file1"
        results_path = "path/to/results"

        expected_output = [
            f"path/to/results/fastqc/unprocessed/file1_1_fastqc.zip",
            f"path/to/results/fastqc/unprocessed/file1_2_fastqc.zip",
            f"path/to/results/fastqc/unprocessed/file1_1_fastqc.html",
            f"path/to/results/fastqc/unprocessed/file1_2_fastqc.html",
            f"path/to/results/fastqc/trimmer/file1_1_fastqc.zip",
            f"path/to/results/fastqc/trimmer/file1_2_fastqc.zip",
            f"path/to/results/fastqc/trimmer/file1_1_fastqc.html",
            f"path/to/results/fastqc/trimmer/file1_2_fastqc.html",
        ]
        assert isinstance(get_fastqc_output(file_name, results_path, mock_pipeline_config), list)
        assert set(expected_output) == (set(get_fastqc_output(file_name, results_path, mock_pipeline_config)))


def test_get_fastqc_output_single_end(mock_pipeline_config):
    mock_pipeline_config.get_config_option.return_value = "trimmer"
    mock_pipeline_config.is_paired_end.return_value = False
    file_name = "file1"
    results_path = "path/to/results"

    expected_output = [
        f"path/to/results/fastqc/unprocessed/file1_fastqc.zip",
        f"path/to/results/fastqc/unprocessed/file1_fastqc.html",
        f"path/to/results/fastqc/trimmer/file1_fastqc.zip",
        f"path/to/results/fastqc/trimmer/file1_fastqc.html",
    ]
    assert isinstance(get_fastqc_output(file_name, results_path, mock_pipeline_config), list)
    assert set(expected_output) == (set(get_fastqc_output(file_name, results_path, mock_pipeline_config)))

