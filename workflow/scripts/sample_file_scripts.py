
import pandas as pd
import pathlib
import json
import os
import re
import yaml

from fetch_data import get_meta_data, get_sra_accessions

#JSON_LOCATION = "/tmp/samples.json"

class SampleFileScripts:
    sample_sheet = None
    paired_end = None
    def __init__(self, config):
        if config["sample_sheet"] is None or config["sample_sheet"] == '':
            self.sample_sheet = "config/samples.csv"
        else:
            self.sample_sheet = config["sample_sheet"]
        self.paired_end = config["paired_end"]

    def make_sample_info(self, json_path) -> dict[str:dict]:
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

        json_dir = "/".join(json_path.split("/")[:-1]) + "/"
        if not os.path.exists(json_dir): os.makedirs(json_dir)
        with open(json_path, "w") as outfile:
            outfile.write(json.dumps(sample_info, indent=4))
        return sample_info

    def is_paired_end(self) -> bool:
        return str(self.paired_end).lower() == "true"

    @staticmethod
    def get_file_info(json_path):
        with open(json_path, "r") as file:
            json_data = json.load(file)
        return json_data

