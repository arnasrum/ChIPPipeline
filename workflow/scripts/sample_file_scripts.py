
import pandas as pd
import pathlib
import json
import os
import re
import yaml

from fetch_data import get_meta_data, get_sra_accessions


JSON_LOCATION="samples.json"

class SampleFileScripts:
    sample_sheet = None
    def __init__(self, config):
        if config["sample_sheet"] is None or config["sample_sheet"] == '':
            self.sample_sheet = "config/samples.csv"
        else:
            self.sample_sheet = config["sample_sheet"]

    def make_sample_info(self) -> dict[str:dict]:
        geo_accession_pattern = re.compile(r"^GSM[0-9]*$")
        geo_accessions = set()
        sample_info: dict = {"public": {}, "provided": {}}
        print(f"sample_sheet: {self.sample_sheet}")
        with open(self.sample_sheet, "r") as file:
            sample_sheet = pd.read_csv(file, keep_default_na=False)
        for index, row in sample_sheet.iterrows():
            sample = row.values[5]
            if geo_accession_pattern.match(sample):
                availability_type = "public"
                geo_accessions.add(sample)
                sample_info[availability_type][sample] = {}
            else:
                # Check if the included files exist
                availability_type = "provided"
                paths = sample.split(";")
                sample = f"{row['mark']}_{row['sample']}_{row['type']}_rep{row['replicate']}".lstrip("_")
                sample_info[availability_type][sample] = {}
                for i, read in enumerate(paths):
                    path = pathlib.Path(read)
                    sample_info[availability_type][sample][f"read{i + 1}"] = {
                        "path": read,
                        "file_extension": "".join(path.suffixes),
                        "file_name": path.name.split("".join(path.suffixes))[0]#.split(f"_{i + 1}")[0]
                    }
                    sample_info[availability_type][sample]['file_name'] = f"{row['mark']}_{row['sample']}_{row['type']}_rep{row['replicate']}".lstrip("_")
            sample_info[availability_type][sample].update({
                "type": row["type"],
                "sample": row["sample"],
                "replicate": row["replicate"],
                "mark": row["mark"],
                "peak_type": row["peak_type"],
            })
        # Handle publicly available files
        fetched_info = get_meta_data(get_sra_accessions(geo_accessions).values())
        sample_info["public"] = {key: value for key, value in map(lambda key: (key, sample_info["public"][key] | fetched_info[key]), sample_info["public"].keys())}

        with open(JSON_LOCATION, "w") as outfile:
            outfile.write(json.dumps(sample_info, indent=4))
        return sample_info


    def generate_json(self):
        pattern = re.compile(r"^GSM[0-9]*$")
        geo_accessions = set()
        sample_info: dict = {}
        print(f"sample_sheet: {self.sample_sheet}")
        with open(self.sample_sheet, "r") as file:
            sample_sheet = pd.read_csv(file, keep_default_na=False)
        for index, row in sample_sheet.iterrows():
            if row["sample"] and row["mark"]:
                sample_info[f"{row['mark']}_{row['sample']}"] = {}
        #for index, row in sample_sheet.iterrows():
            #for sample in sample_info.keys():
                #if row["mark"] in sample:
        print(sample_info)

    @staticmethod
    def get_file_info():
        with open(JSON_LOCATION, "r") as file:
            json_data = json.load(file)
        return json_data

    @staticmethod
    def get_macs_input() -> dict[str: dict]:
        with open(JSON_LOCATION) as file:
            samples = json.load(file)
        macs_input = {}
        control_files: list[(str, int, str)] = []
        for key, entry in SampleFileScripts.__flatten_dict(samples).items():
            replicate = str(entry["replicate"])
            if entry["type"] == "treatment":
                file_name = f"{entry['mark']}_{entry['sample']}"
                if not file_name in macs_input: macs_input[file_name] = {}
                if not replicate in macs_input[file_name]:
                    macs_input[file_name][replicate] = {"treatment": [entry['file_name']], "control": [], "peak_type": entry['peak_type']}
                else:
                    macs_input[file_name][replicate]["treatment"].append(entry['file_name'])
            elif entry["type"] == "control":
                control_files.append((entry['sample'], replicate, entry['file_name']))
            else:
                raise Exception(f"Entry type unrecognized for {entry['file_name']}")
        for control in control_files:
            for file_name in filter(lambda file_name: control[0] in file_name, macs_input.keys()):
                macs_input[file_name][control[1]]["control"].append(control[2])

        invalid_input: list[(str, str)] = []
        for file_name in macs_input:
            for missing_pair_replicate in filter(lambda replicate: len(macs_input[file_name][replicate]["treatment"]) == 0 or len(macs_input[file_name][replicate]["control"]) == 0, macs_input[file_name].keys()):
                invalid_input.append((file_name, missing_pair_replicate))
        for file, replicate in invalid_input:
            macs_input[file].pop(replicate)
        return macs_input

    @staticmethod
    def __flatten_dict(old_dict: dict) -> dict:
        new_dict = {}
        for key, value in old_dict.items():
            new_dict = new_dict | value
        return new_dict

def is_paired_end() -> bool:
    with open ("config/config.yml") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    if str(config["paired_end"]).lower() == "true":
        return True
    else:
        return False


if __name__ == "__main__":
    config = yaml.load(open("config/config.yml", "r"), Loader=yaml.FullLoader)
    config["sample_sheet"] = "config/samplesTest.csv"
    SampleFileScripts(config).make_sample_info()
    print(SampleFileScripts(config).get_macs_input())
