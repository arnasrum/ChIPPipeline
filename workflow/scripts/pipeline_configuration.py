import os
import re
import json
import pathlib
import pandas as pd
import numpy as np
from workflow.scripts.fetch_data import get_meta_data, get_sra_accessions

class InputException(Exception):
    """Custom exception class for handling input-related errors."""
    pass

class PipelineConfiguration:
    """
    Handles the configuration and sample information for the pipeline.

    This class is responsible for parsing the configuration yaml file, validating
    the sample sheet and configuration arguments, and generating a dictionary of
    parsed sample data used throughout the pipeline.
    """
    def __init__(self, config: dict):
        """
        Initializes the PipelineConfiguration with the provided configuration dictionary.

        Args:
            config: A dictionary containing the pipeline configuration.
        """
        self.sample_info = None
        if config["sample_sheet"] is None or config["sample_sheet"] == '':
            self.sample_sheet = "config/samples.csv"
        else:
            self.sample_sheet = config["sample_sheet"]
        self.config = config

    def __validate_config(self):
        """
        Validates the main configuration parameters and arguments.

        Checks if specified boolean arguments have valid 'true' or 'false' values.
        Also validates potentially injectable arguments within the configuration
        dictionary for invalid characters to prevent command injection.

        Raises:
            InputException: If a boolean configuration argument has an invalid value
                            or if any configuration argument contains invalid characters.
        """
        boolean_keys = ["generate_fastqc_reports"]
        for key in boolean_keys:
            if not (str(self.config[key]).lower() == "true" or str(self.config[key]).lower() == "false"):
                raise InputException(f"The configuration argument; {key}, must be true or false.")
        char_whitelist = re.compile(r"[a-zA-Z0-9,._:\s\"\'\-]")
        injectable_options = ["args", "run_options", "mode"]
        for tool in self.config:
            for key in injectable_options:
                if isinstance(self.config[tool], dict) and key in self.config[tool].keys() and self.config[tool][key]:
                    if next(filter(lambda char: not re.match(char_whitelist, char), self.config[tool][key]), None):
                        raise InputException(f"The configuration argument; {tool} {key}, contains invalid characters.")



    def make_sample_info(self) -> dict[str:dict]:
        """
        Creates a dictionary of with parsed sample information from the sample sheet.

        Validates the configuration and sample sheet, then processes the sample sheet
        to create a nested dictionary structure separating public and provided data.
        Fetches metadata for public accessions and saves the combined sample info.

        Returns:
            A dictionary containing the structured sample information.
        """

        geo_accession_pattern = re.compile(r"^GSM[0-9]*$")
        geo_accessions = set()
        sample_info: dict = {"public": {}, "provided": {}}

        if not os.path.isfile(self.sample_sheet):
            raise InputException(f"Sample sheet; \"{self.sample_sheet}\" does not exist.")
        if self.sample_sheet.endswith(".csv"):
            with open(self.sample_sheet, "r") as file:
                sample_sheet = pd.read_csv(file, keep_default_na=False)
        elif self.sample_sheet.endswith(".json"):
            with open(self.sample_sheet, "r") as file:
                sample_sheet = pd.read_json(file)
        else:
            raise Exception("Sample sheet file format not supported. Only .csv and .json are supported.")
        sample_sheet = sample_sheet.fillna(np.nan).replace([np.nan], [None])
        PipelineConfiguration.__validate_sample_sheet(sample_sheet)
        self.__validate_config()

        for index, row in sample_sheet.iterrows():
            sample = ""
            availability_type = ""
            if (row["file_path"] is None or row["file_path"] == "") and geo_accession_pattern.match(row["accession"]):
                availability_type = "public"
                sample = row["accession"]
                geo_accessions.add(row["accession"])
                sample_info[availability_type][sample] = {}
            if row["file_path"]:
                availability_type = "provided"
                paths = row["file_path"].split(";")
                genome_name = row['genome'].split('/')[-1].split('.')[0]
                sample = f"{row['sample']}_{row['type']}_rep{row['replicate']}_{genome_name}".lstrip('_')
                if row["mark"]:
                    sample = f"{row['mark']}_{sample}"
                if row["accession"]:
                    sample = f"{row['accession']}_{sample}"
                sample_info[availability_type][sample] = {}
                for i, read in enumerate(paths):
                    path = pathlib.Path(read)
                    sample_info[availability_type][sample][f"read{i + 1}"] = {
                        "path": read,
                        "file_extension": "".join(path.suffixes),
                        "file_name": path.name.split("".join(path.suffixes))[0]
                    }
                    sample_info[availability_type][sample]['file_name'] = sample
            sample_info[availability_type][sample].update({
                "type": row["type"],
                "sample": row["sample"],
                "replicate": row["replicate"],
                "mark": row["mark"],
                "peak_type": row["peak_type"],
                "genome": row["genome"],
                "paired_end": str(row["paired_end"]).lower(),
            })
        # Handle publicly available files
        fetched_info = get_meta_data(get_sra_accessions(geo_accessions).values())
        sample_info["public"] = {key: value for key, value in map(lambda key: (key, sample_info["public"][key] | fetched_info[key]), sample_info["public"].keys())}
        for sample in sample_info['public']:
            sample_info['public'][sample]['file_name'] += f"_{sample_info['public'][sample]['genome'].split('/')[-1].split('.')[0]}"

        self.sample_info = sample_info
        return sample_info

    def is_paired_end(self, file_name: str) -> bool:
        """
        Checks if a sample is paired-end based on its file name.

        Args:
            file_name: The file name of the sample to check.

        Returns:
            True if the sample is paired-end, False otherwise.
        """
        return str(self.get_sample_entry_by_file_name(file_name)["paired_end"]).lower() in ["true", "yes", "y", "t"]

    def get_all_file_names(self) -> list[str]:
        """
        Retrieves all generated file names from the sample information.

        Returns:
            A list of all file names for the samples.
        """
        return [sample_info["file_name"] for sample_info in PipelineConfiguration.__flatten_dict(self.sample_info).values()]

    def __group_control_files(self) -> dict[str, dict[str, list[str]]]:
        """
        Groups control files by sample and replicate.

        Returns:
            A nested dictionary where the first key is the sample group
            (sample_name_genome) and the second key is the replicate,
            containing a list of control file names.
        """
        control_files: dict[str, dict[str, list[str]]] = {}
        for key, entry in self.__flatten_dict(self.sample_info).items():
            if entry['type'] != "control": continue
            group_name = f"{entry['sample']}_{self.get_sample_genome(entry['file_name'])}"
            if not group_name in control_files:
                control_files[group_name] = {}
            if not entry['replicate'] in control_files[group_name]:
                control_files[group_name][entry['replicate']] = []
            control_files[group_name][entry['replicate']].append(entry['file_name'])
        return control_files

    def get_control_files(self, treatment_file: str) -> list[str]:
        """
        Fetches the control files corresponding to a given treatment file.

        Finds the control samples that match the sample name, genome, and
        replicate of the provided treatment file.

        Args:
            treatment_file: The file name of the treatment sample.

        Returns:
            A list of control file names associated with the treatment sample.

        Raises:
            Exception: If the provided treatment file does not exist in the sample info.
        """
        treatment_entry = next(filter(lambda entry: entry["file_name"] == treatment_file, self.__flatten_dict(self.sample_info).values()))
        if treatment_entry is None:
            raise Exception(f"File: {treatment_file}; does not exists in sample info")
        sample = f"{treatment_entry['sample']}_{self.get_sample_genome(treatment_file)}"
        replicate = treatment_entry["replicate"]
        control_files = self.__group_control_files().get(sample, {}).get(replicate, [])
        return control_files

    def group_treatment_files(self) -> dict[str, dict[str, list[str]]]:
        """
        Groups treatment files by mark, sample, and genome.

        Returns:
            A nested dictionary where the first key is the treatment group
            (mark_sample_genome) and the second key is the replicate,
            containing a list of treatment file names.
        """
        treatment_samples: dict[str, dict[str, list[str]]] = {}
        for key, entry in self.__flatten_dict(self.sample_info).items():
            if entry['type'] != "treatment": continue
            genome = entry["genome"].split("/")[-1].split(".")[0]
            group = f"{entry['mark']}_{entry['sample']}_{genome}"
            if not group in treatment_samples:
                treatment_samples[group] = {}
            if not entry['replicate'] in treatment_samples[group]:
                treatment_samples[group][entry['replicate']] = []
            treatment_samples[group][entry['replicate']].append(entry['file_name'])
        return treatment_samples

    def get_treatment_files(self, return_entries = False) -> list[str]:
        """
        Retrieves a list of treatment files.

        Args:
            return_entries: If True, returns the full sample information
                            dictionaries for treatment samples. If False
                            (default), returns only the file names.

        Returns:
            A list of treatment file names or sample information dictionaries.
        """
        if return_entries:
            return [*filter(lambda entry: entry["type"] == "treatment", self.__flatten_dict(self.sample_info).values())]
        else:
            return [entry['file_name'] for entry in filter(lambda entry: entry["type"] == "treatment", self.__flatten_dict(self.sample_info).values())]

    def get_sample_entry_by_file_name(self, file_name: str) -> dict:
        """
        Retrieves the dictionary entry for a sample based on the given file name.

        Args:
            file_name: The file name of the sample.

        Returns:
            The dictionary containing all information for the specified sample.
        """
        return next(filter(lambda entry: entry["file_name"] == file_name, self.__flatten_dict(self.sample_info).values()))

    def get_sample_genome(self, file_name: str) -> str:
        """
        Fetches the genome associated with a sample based on the given file name.

        In case the genome was provided as a file path, the function extracts
        the file name of the genome.

        Args:
            file_name: The file name of the sample.

        Returns:
            The extracted genome name.
        """
        genome = self.get_sample_entry_by_file_name(file_name)['genome'].split("/")[-1].split(".")[0]
        return genome

    @staticmethod
    def __validate_sample_sheet(sample_sheet: pd.DataFrame):
        for index, row in sample_sheet.iterrows():
            if not row["type"] in ["treatment", "control"]:
                raise InputException(f"Row {index} in \"type\" column contains invalid value. {row['type']}")
            if  row["type"] == "treatment" and not row["peak_type"] in ["narrow", "broad", ""]:
                if row["peak_type"] == "" or row["peak_type"] is None:
                    raise InputException(f"Row {index} in \"peak_type\" treatment samples must have defined peak type.")
                raise InputException(f"Row {index} in \"peak_type\" column contains invalid value. {row['peak_type']}")
            if not isinstance(row["replicate"], int) or int(row["replicate"]) < 0:
                # something wrong here, if genome is defined this gives the wrong message
                raise InputException(f"Row {index} in \"replicate\" column contains invalid value. {row['replicate']} must be a positive integer.")
            if not row["accession"] and not row["file_path"]:
                raise InputException(f"Row {index} in column \"accession\" contains invalid value.")
            if row["file_path"]:
                for file in row["file_path"].split(";"):
                    if not os.path.isfile(file):
                        raise InputException(f"Provided file: {file}, in row {index}, does not exist.")
            if row["file_path"] and str(row["paired_end"]).lower() == 'true' and len(row["file_path"].split(";")) == 1:
                raise InputException(f"Running pipeline in paired end mode, but only one read provided for row {index} in sample sheet.")

    @staticmethod
    def __flatten_dict(old_dict: dict) -> dict:
        """
        Flattens a nested dictionary (specifically the sample_info structure).

        Combines the 'public' and 'provided' nested dictionaries into a single
        level dictionary.

        Args:
            old_dict: The nested dictionary to flatten (e.g., self.sample_info).

        Returns:
            A flattened dictionary.
        """
        new_dict = {}
        for key, value in old_dict.items():
            new_dict = new_dict | value
        return new_dict
