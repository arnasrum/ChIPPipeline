import pandas as pd
import json

def check_diff() -> bool:
    samples_csv = open("config/samples.csv")
    json_file = open("config/samples.json")
    sample_sheet = pd.read_csv(samples_csv, keep_default_na=False)
    samples_json = json.load(json_file)
    samples_csv.close()
    json_file.close()
    if len(sample_sheet) != len(samples_json["public"]) + len(samples_json["provided"]): return False
    for index, row in sample_sheet.iterrows():
        if row["accession"] in samples_json["public"]:
            print(row["accession"], "in public")
            for header in sample_sheet.head():
                if header == "accession": continue
                if not header in samples_json["public"][row["accession"]]: return False
            continue
        if row["accession"] in samples_json["provided"]:
            print(row["accession"], "in provided")
            continue
        return False
    return True
