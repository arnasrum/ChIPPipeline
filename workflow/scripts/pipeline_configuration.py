import os
import re
import json
import yaml
import pathlib
import pandas as pd
from workflow.scripts.fetch_data import get_meta_data, get_sra_accessions

class InputException(Exception):
    pass

class PipelineConfiguration:
    def __init__(self, config: dict):
        self.sample_info = None
        if config["sample_sheet"] is None or config["sample_sheet"] == '':
            self.sample_sheet = "config/samples.csv"
        else:
            self.sample_sheet = config["sample_sheet"]
        self.paired_end = config["paired_end"]
        self.json_path = config["json_path"]
        self.config = config

    def __verify_config(self):
        config = self.config
        if config["genome"] == "" or config["genome"] is None:
            raise InputException("Please specify a genome in the configuration file.")
        # Verify boolean arguments
        boolean_keys = ["paired_end", "generate_fastqc_reports"]
        for key in boolean_keys:
            if not (str(config[key]).lower() == "true" or str(config[key]).lower() == "false"):
                raise InputException(f"The configuration argument; {key}, must be true or false.")

    def __verify_sample_sheet(self):
        with open(self.sample_sheet, "r") as file:
            sample_sheet = pd.read_csv(file, keep_default_na=False)
        for index, row in sample_sheet.iterrows():
            if not row["type"] in ["treatment", "control"]:
                raise InputException(f"Row {index} in \"type\" column contains invalid value. {row['type']}")
            if not row["peak_type"] in ["narrow", "broad", ""]:
                if row["type"] == "treatment" and row["peak_type"] == "":
                    raise InputException(f"Row {index} in \"peak_type\" treatment samples must have defined peak type.")
                raise InputException(f"Row {index} in \"peak_type\" column contains invalid value. {row['type']}")
            if not isinstance(row["replicate"], int):
                # something wrong here, if genome is defined this gives the wrong message
                raise InputException(f"Row {index} in \"replicate\" column contains invalid value. {row['replicate']} must be an integer.")
            if not row["accession"]:
                raise InputException(f"Row {index} in column \"accession\" contains invalid value.")
            if row["file_path"]:
                for file in row["file_path"].split(";"):
                    if not os.path.isfile(file):
                        raise InputException(f"Provided file: {file}, in row {index}, does not exist.")
            if row["file_path"] and self.is_paired_end() and len(row["file_path"].split(";")) == 1:
                raise InputException(f"Running pipeline in paired end mode, but only one read provided for row {index} in sample sheet.")

    def make_sample_info(self) -> dict[str:dict]:
        self.__verify_sample_sheet()
        self.__verify_config()
        geo_accession_pattern = re.compile(r"^GSM[0-9]*$")
        geo_accessions = set()
        sample_info: dict = {"public": {}, "provided": {}}
        #print(f"sample_sheet: {self.sample_sheet}")
        with open(self.sample_sheet, "r") as file:
            sample_sheet = pd.read_csv(file, keep_default_na=False)
        for index, row in sample_sheet.iterrows():
            sample = ""
            availability_type = ""
            if row["file_path"] == "" and geo_accession_pattern.match(row["accession"]):
                availability_type = "public"
                sample = row["accession"]
                geo_accessions.add(row["accession"])
                sample_info[availability_type][sample] = {}
            if row["file_path"]:
                availability_type = "provided"
                paths = row["file_path"].split(";")
                sample = f"{row['mark']}_{row['sample']}_{row['type']}_rep{row['replicate']}_{row['genome'].split('/')[-1].split('.')[0]}"
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
            })
        # Handle publicly available files
        fetched_info = get_meta_data(get_sra_accessions(geo_accessions).values())
        sample_info["public"] = {key: value for key, value in map(lambda key: (key, sample_info["public"][key] | fetched_info[key]), sample_info["public"].keys())}
        for sample in sample_info['public']:
            sample_info['public'][sample]['file_name'] += f"_{sample_info['public'][sample]['genome'].split('/')[-1].split('.')[0]}"

        json_dir = "/".join(self.json_path.split("/")[:-1]) + "/"
        if not os.path.exists(json_dir): os.makedirs(json_dir)
        with open(self.json_path, "w") as outfile:
            outfile.write(json.dumps(sample_info, indent=4))
        self.sample_info = sample_info
        return sample_info

    def is_paired_end(self) -> bool:
        return str(self.paired_end).lower() == "true"

    def check_diff(self) -> bool:
        samples_csv = open(self.sample_sheet)
        json_file = open(self.json_path)
        sample_sheet = pd.read_csv(samples_csv, keep_default_na=False); samples_csv.close()
        samples_json = json.load(json_file); json_file.close()
        samples_json = PipelineConfiguration.__flatten_dict(samples_json)
        if len(sample_sheet) != len(samples_json["public"]) + len(samples_json["provided"]): return False
        for index, row in sample_sheet.iterrows():
            if row["accession"] in samples_json["public"]:
                # print(row["accession"], "in public")
                for header in sample_sheet.head():
                    if header == "accession": continue
                    if not header in samples_json["public"][row["accession"]]: return False
                    if not row[header] == samples_json["public"][row["accession"]][header]: return False
                continue
            if row["accession"] in samples_json["provided"]:

                # print(row["accession"], "in provided")
                continue
            return False
        return True

    def get_all_file_names(self) -> list[str]:
        return [sample_info["file_name"] for sample_info in PipelineConfiguration.__flatten_dict(self.sample_info).values()]


    def __group_control_files(self) -> dict[str, dict[str, list[str]]]:
        control_files: dict[str, dict[str, list[str]]] = {}
        for key, entry in self.__flatten_dict(self.sample_info).items():
            if entry['type'] != "control": continue
            if not entry['sample'] in control_files:
                control_files[entry['sample']] = {}
            if not entry['replicate'] in control_files[entry['sample']]:
                control_files[entry['sample']][entry['replicate']] = []
            control_files[entry['sample']][entry['replicate']].append(entry['file_name'])
        return control_files

    def get_control_files(self, treatment_file: str) -> list[str]:
        treatment_entry = next(filter(lambda entry: entry["file_name"] == treatment_file, self.__flatten_dict(self.sample_info).values()))
        if treatment_entry is None:
            raise Exception(f"File: {treatment_file}; does not exists in sample info")
        sample = treatment_entry["sample"]
        replicate = treatment_entry["replicate"]
        control_files = self.__group_control_files().get(sample, {}).get(replicate, [])
        return control_files

    def group_treatment_files(self) -> dict[str, dict[str, list[str]]]:
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
        if return_entries:
            return [*filter(lambda entry: entry["type"] == "treatment", self.__flatten_dict(self.sample_info).values())]
        else:
            return [entry['file_name'] for entry in filter(lambda entry: entry["type"] == "treatment", self.__flatten_dict(self.sample_info).values())]

    def get_genome_code(self):
        return self.config["genome"].split("/")[-1].split(".")[0]

    def get_sample_entry_by_file_name(self, file_name: str) -> dict:
        return next(filter(lambda entry: entry["file_name"] == file_name, self.__flatten_dict(self.sample_info).values()))

    def get_sample_genome(self, file_name: str) -> str:
        genome = self.get_sample_entry_by_file_name(file_name)['genome'].split("/")[-1].split(".")[0]
        return genome

    @staticmethod
    def __flatten_dict(old_dict: dict) -> dict:
        new_dict = {}
        for key, value in old_dict.items():
            new_dict = new_dict | value
        return new_dict


if __name__ == "__main__":
    file = open("config/config.yml")
    config = yaml.load(file, Loader=yaml.FullLoader)
    file.close()
    sfs = PipelineConfiguration(config)
    sfs.make_sample_info()
    #sfs.group_samples()
    sfs.get_control_files("GSM1871972_H3K4me3_8cell_treatment_rep1")
