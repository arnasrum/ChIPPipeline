from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory, _TemporaryFileWrapper
from typing import Iterator, List
import pandas as pd
from workflow.scripts.pipeline_configuration import PipelineConfiguration
import pytest
import json

@pytest.fixture
def mock_sample_sheet_json() -> Iterator[NamedTemporaryFile]:
    pre_created_temp_files: List[NamedTemporaryFile] = []
    temp_read1 = NamedTemporaryFile(mode="w+", suffix=".fq.gz")
    temp_read1.flush() # Ensure the empty file is created on disk
    pre_created_temp_files.append(temp_read1)
    temp_read2 = NamedTemporaryFile(mode="w+", suffix=".fq.gz")
    temp_read2.flush()
    pre_created_temp_files.append(temp_read2)

    sample_sheet_data = [
        {
            "mark": "antibody",
            "sample": "sample1",
            "type": "treatment",
            "replicate": 1,
            "peak_type": "narrow",
            "paired_end": True,
            "accession": None,
            "file_path": f"{temp_read1.name};{temp_read2.name}",
            "genome": "data/genome.fa.gz"
        },
        {
            "mark": "antibody",
            "sample": "sample1",
            "type": "treatment",
            "replicate": 1,
            "peak_type": "narrow",
            "paired_end": "true",
            "accession": "GSM123",
            "file_path": f"{temp_read1.name};{temp_read2.name}",
            "genome": "data/genome.fa.gz"
        },
        {
            "mark": None,
            "sample": "sample1",
            "type": "control",
            "replicate": 1,
            "peak_type": None,
            "paired_end": False,
            "accession": None,
            "file_path": f"{temp_read1.name};{temp_read2.name}",
            "genome": "data/genome.fa.gz"
        },
        {
            "mark": "antibody",
            "sample": "sample1",
            "type": "treatment",
            "replicate": "2",
            "peak_type": "narrow",
            "paired_end": False,
            "accession": "GSM123451",
            "file_path": temp_read1.name,
            "genome": "genome"
        }
    ]

    with NamedTemporaryFile(suffix='.json', mode="w+") as file:
        json.dump(sample_sheet_data, file, indent=4)
        file.flush()
        yield file


@pytest.fixture
def pipeline_config(mock_sample_sheet_json) -> PipelineConfiguration:
    config = {"sample_sheet": mock_sample_sheet_json.name, "generate_fastqc_reports": True}
    pipeline_config = PipelineConfiguration(config)
    pipeline_config.parse_sample_sheet()
    return pipeline_config

def test_pipeline_configuration_parse_sample_sheet(pipeline_config):
    assert "public" in pipeline_config.sample_info and "provided" in pipeline_config.sample_info
    assert "antibody_sample1_rep1_genome_treatment" in pipeline_config.sample_info["provided"] and \
        "antibody_sample1_rep2_genome_treatment" in pipeline_config.sample_info["provided"] and \
        "sample1_rep1_genome_control" in pipeline_config.sample_info["provided"] and \
        len(pipeline_config.sample_info["provided"]["antibody_sample1_rep1_genome_treatment"]) == 2

    entry = pipeline_config.sample_info["provided"]["antibody_sample1_rep1_genome_treatment"][0]
    assert "read1" in entry and "read2" in entry and "file_name" in entry  and"sample" in entry and \
            "replicate" in entry and "genome" in entry and "paired_end" in entry and "peak_type" in entry and \
            entry["peak_type"] == "narrow" and entry["replicate"] == 1 and str(entry["paired_end"]).lower() == "true"

    entry = pipeline_config.sample_info["provided"]["antibody_sample1_rep2_genome_treatment"][0]
    assert "read1" in entry and not "read2" in entry


def test_pipeline_configuration_get_all_file_names(pipeline_config):
    expected_output = {"antibody_sample1_treatment_rep1_genome",
                       "GSM123_antibody_sample1_treatment_rep1_genome_pooled2",
                       "GSM123451_antibody_sample1_treatment_rep2_genome",
                       "sample1_control_rep1_genome"
    }
    assert expected_output == set(pipeline_config.get_all_file_names())

def test_pipeline_configuration_get_treatment_files(pipeline_config):
    # No pooled
    expected_output = {"antibody_sample1_treatment_rep1_genome",
                       "GSM123451_antibody_sample1_treatment_rep2_genome",
    }
    assert expected_output == set(pipeline_config.get_treatment_files())



def test_pipeline_configuration_is_paired_end(pipeline_config):
    assert pipeline_config.is_paired_end("antibody_sample1_treatment_rep1_genome") == True
    assert pipeline_config.is_paired_end("GSM123_antibody_sample1_treatment_rep1_genome_pooled2") == True
    assert pipeline_config.is_paired_end("sample1_control_rep1_genome") == False

def test_pipeline_configuration_get_sample_entry_by_file_name(pipeline_config):
    groups = pipeline_config.group_treatment_files()
    assert "antibody_sample1_genome" in groups and 1 in groups["antibody_sample1_genome"] and 2 in groups["antibody_sample1_genome"] and \
        "antibody_sample1_treatment_rep1_genome" == groups["antibody_sample1_genome"][1][0] and \
        "GSM123_antibody_sample1_treatment_rep1_genome_pooled2" == groups["antibody_sample1_genome"][1][1] and \
        "GSM123451_antibody_sample1_treatment_rep2_genome" == groups["antibody_sample1_genome"][2][0]

