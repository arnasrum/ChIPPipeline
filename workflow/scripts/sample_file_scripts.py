
import pandas as pd
import pathlib
import json
import os
import re
import yaml

from fetch_data import get_meta_data, get_sra_accessions
from input_scripts import symlink_input


class SampleFileScripts:
    def __init__(self, config: dict):
        if config["sample_sheet"] is None or config["sample_sheet"] == '':
            self.sample_sheet = "config/samples.csv"
        else:
            self.sample_sheet = config["sample_sheet"]
        self.paired_end = config["paired_end"]
        self.json_path = config["json_path"]

    def make_sample_info(self) -> dict[str:dict]:
        geo_accession_pattern = re.compile(r"^GSM[0-9]*$")
        geo_accessions = set()
        sample_info: dict = {"public": {}, "provided": {}}
        print(f"sample_sheet: {self.sample_sheet}")
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
                # Check if the included files exist
                availability_type = "provided"
                paths = row["file_path"].split(";")
                sample = f"{row['mark']}_{row['sample']}_{row['type']}_rep{row['replicate']}".lstrip("_")
                if row["accession"]:
                    sample = f"{row['accession']}_{sample}"
                sample_info[availability_type][sample] = {}
                for i, read in enumerate(paths):
                    path = pathlib.Path(read)
                    sample_info[availability_type][sample][f"read{i + 1}"] = {
                        "path": read,
                        "file_extension": "".join(path.suffixes),
                        "file_name": path.name.split("".join(path.suffixes))[0]#.split(f"_{i + 1}")[0]
                    }
                    sample_info[availability_type][sample]['file_name'] = sample
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

        json_dir = "/".join(self.json_path.split("/")[:-1]) + "/"
        if not os.path.exists(json_dir): os.makedirs(json_dir)
        with open(self.json_path, "w") as outfile:
            outfile.write(json.dumps(sample_info, indent=4))
        return sample_info

    def is_paired_end(self) -> bool:
        return str(self.paired_end).lower() == "true"

    def check_diff(self) -> bool:
        samples_csv = open(self.sample_sheet)
        json_file = open(self.json_path)
        sample_sheet = pd.read_csv(samples_csv, keep_default_na=False); samples_csv.close()
        samples_json = json.load(json_file); json_file.close()
        samples_json = SampleFileScripts.__flatten_dict(samples_json)
        print(samples_json)
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

    @staticmethod
    def get_file_info(json_path):
        with open(json_path, "r") as file:
            json_data = json.load(file)
        return json_data

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
    sfs = SampleFileScripts(config)
    sfs.make_sample_info("")
    symlink_input(config["json_path"], "test")



