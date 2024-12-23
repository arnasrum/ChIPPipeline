
import pandas as pd
import pathlib
import json
import os
import re
import yaml

from fetch_data import get_meta_data, get_sra_accessions

def find_prefix(strings):
    if not strings:
        return ""
    if len(strings) == 1:
        return strings[0]
    i = 0; j = 0
    prefix = [] 
    while len(strings[0]) > i and len(strings[1]) > j and strings[0][i] == strings[1][j]:
        prefix.append(strings[0][i])
        i+=1; j+=1
    return "".join(prefix)[:-1]

class SampleFileScripts:

    sample_sheet = None
    def __init__(self, config):
        if config["sample_sheet"] is None or config["sample_sheet"] == '':
            self.sample_sheet = "config/samples.csv"
        else:
            self.sample_sheet = config["sample_sheet"]

    def make_sample_info(self) -> dict[str:dict]:
        pattern = re.compile(r"^GSM[0-9]*$")
        geo_accessions = set()
        sample_info: dict = {"public": {}, "provided": {}}
        with open(self.sample_sheet, "r") as file:
            for index, row in pd.read_csv(file, keep_default_na=False).iterrows():
                sample = row.values[5]
                if pattern.match(sample):
                    availability_type = "public"
                    geo_accessions.add(sample)
                    sample_info[availability_type][sample] = {}
                else:
                    # Check if the included files exist
                    reads = sample.split(";")
                    availability_type = "provided"
                    sample_info[availability_type][sample] = {}
                    read_file_names = []
                    for i, read in enumerate(reads):
                        path = pathlib.Path(read)
                        sample_info[availability_type][sample][f"read{i + 1}"] = {
                            "path": read,
                            "file_extension": "".join(path.suffixes),
                            "file_name": path.name.split("".join(path.suffixes))[0]#.split(f"_{i + 1}")[0]
                        }
                        read_file_names.append(sample_info[availability_type][sample][f"read{i + 1}"]["file_name"])
                    sample_info[availability_type][sample]["cleanFileName"] = find_prefix(read_file_names)
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

        with open(f"config/samples.json", "w") as outfile:
            outfile.write(json.dumps(sample_info, indent=4))
        return sample_info

    @staticmethod
    def get_file_info():
        with open("config/samples.json", "r") as file:
            json_data = json.load(file)
        return json_data

    @staticmethod
    def get_macs_input() -> dict[str: dict]:
        with open("config/samples.json") as file: samples = json.load(file)
        macs_input = {}
        control_files: list[(str, int, str)] = []
        for entry in SampleFileScripts.__flatten_dict(samples).values():
            replicate = str(entry["replicate"])
            if entry["type"] == "treatment":
                file_name = f"{entry['mark']}_{entry['sample']}"
                if not file_name in macs_input: macs_input[file_name] = {}
                if not replicate in macs_input[file_name]:
                    macs_input[file_name][replicate] = {"treatment": [entry['cleanFileName']], "control": [], "peak_type": entry["peak_type"]}
                else:
                    macs_input[file_name][replicate]["treatment"].append(entry['cleanFileName'])
            elif entry["type"] == "control":
                control_files.append((entry['sample'], replicate, entry['cleanFileName']))
            else:
                raise Exception(f"Entry type unrecognized for {entry['cleanFileName']}")
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

if __name__ == "__main__":
    config = yaml.load(open("config/config.yml", "r"), Loader=yaml.FullLoader)
    SampleFileScripts(config).make_sample_info()