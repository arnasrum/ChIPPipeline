from requests import get, Response
from xml.etree import ElementTree
import pandas as pd
import pathlib
import time
import json
import os
import re

REQUEST_TIMEOUT: int = 30
NUM_FAILED_REQUESTS_ALLOWED: int = 5

def make_sample_info(sample_sheet:str= "config/samples.csv") -> dict[str:dict]:

    #if os.path.isfile("config/samples.json"):
        #with open("config/samples.json", "r") as file:
            #return json.load(file)

    pattern = re.compile(r"^GSM[0-9]*$")
    geo_accessions = set()
    sample_info: dict = {"public": {}, "provided": {}}
    with open(sample_sheet, "r") as file:
        for index, row in pd.read_csv(file).iterrows():
            for sample in row.values[5:]:
                if pattern.match(sample):
                    availability_type = "public"
                    geo_accessions.add(sample)
                    sample_info[availability_type][sample] = {}
                else:
                    availability_type = "provided"
                    path = pathlib.Path(sample)
                    file_extension = "".join(path.suffixes)
                    file_name = path.name.split(file_extension)[0]
                    sample_info[availability_type][file_name] = {}
                    sample_info[availability_type][file_name]["cleanFileName"] = file_name
                    sample_info[availability_type][file_name]["fileExtension"] = file_extension
                    sample_info[availability_type][file_name]["path"] = sample.split(file_name + file_extension)[0]
                    sample = file_name
                sample_info[availability_type][sample]["type"] = row["type"]
                sample_info[availability_type][sample]["sample"] = row["sample"]
                sample_info[availability_type][sample]["replicate"] = row["replicate"]
                sample_info[availability_type][sample]["mark"] = row["mark"]
                sample_info[availability_type][sample]["peak_type"] = row["peak_type"]


    # Handle publicly available files
    fetched_info = get_meta_data(get_sra_accessions(geo_accessions).values())
    sample_info["public"] = {key: value for key, value in map(lambda key: (key, sample_info["public"][key] | fetched_info[key]), sample_info["public"].keys())}

    with open(f"config/samples.json", "w") as outfile:
        outfile.write(json.dumps(sample_info, indent=4))
    return sample_info

def get_all_sample_file_paths(config: dict, include_directories=True) -> list[str]:
    directory = "resources/reads/" if include_directories else ""
    file_paths = []
    with open("config/samples.json", "r") as file:
        data: dict[str:dict] = json.load(file)
        for file_type in data:
            for file_info in data[file_type].values():
                path = f"{directory}{file_info['cleanFileName']}"
                file_paths.append(f"{path}_1.fastq")
                if config["paired_end"]: file_paths.append(f"{path}_2.fastq")

    return file_paths


def get_sra_accessions(geo_accessions: set[str]) -> dict[str:str]:
    '''
        Given a list of GEO accessions retrieve their related SRA accessions 
    '''
    if len(geo_accessions) == 0: return dict()
    enterez_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term="
    enterez_url += "+OR+".join(geo_accessions)
    print(enterez_url)
    response: Response = __poll_request(enterez_url)
    xml_response: ElementTree.Element = ElementTree.fromstring(response.content)
    id_list = list(filter(lambda item: int(item) > 299999999, [sample_id.text for sample_id in xml_response[3]]))
    enterez_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&id=" + ",".join(id_list)
    print(enterez_url)
    #response = get(enterez_url, timeout=REQUEST_TIMEOUT)
    response = __poll_request(enterez_url)
    response_text: str = response.content.decode().lstrip("\n").rstrip("\n")
    gsm_to_sra_map: dict[str:str] = {}
    response_text = response_text.replace("\n\n", "\t;:.,").replace("\n", "")
    for match in response_text.split("\t;:.,"):
        sub_response: str = match
        accession_position: re.Match = re.search(r"Accession:\sGSM[0-9]*", sub_response)
        sra_accession_position: re.Match = re.search(r"SRX[0-9]*", sub_response)
        accession: str = sub_response[accession_position.span()[0]: accession_position.span()[1]]
        sra_accession: str = sub_response[sra_accession_position.span()[0]: sra_accession_position.span()[1]]
        gsm_to_sra_map[accession] = sra_accession
    return gsm_to_sra_map

def get_meta_data(sra_accessions: list[str]) -> dict[str: dict]:
    '''
        Fetch metadata for SRA samples 
    '''
    if len(sra_accessions) == 0:
        return dict()
    enterez_url: str = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={','.join(sra_accessions)}"
    print(enterez_url)
    #response: Response = get(enterez_url, timeout=REQUEST_TIMEOUT)
    response = __poll_request(enterez_url)
    root: ElementTree.Element = ElementTree.fromstring(response.content)
    meta_data: dict[str: dict] = {}
    for node in root:
        run_accessions: list[str] = []
        geo_accession = node[4].attrib["alias"]
        meta_data[geo_accession] = {"cleanFileName": __make_clean_file_name(node[0][1].text)}
        for run in node[6]:
            run_accessions.append(run.attrib["accession"])
        meta_data[geo_accession]["runs"] = run_accessions
    #response = xmltodict.parse(result.content)
    return meta_data

def get_file_names(include_provided: bool = True, include_publicly_available: bool = True) -> list[str]:
    """
        Returns a list of file names 
    """
    if not include_provided and not include_publicly_available:
        return list()

    if not os.path.isfile("config/samples.json"):
        make_sample_info()
    #makeSampleInfo()
    with open("config/samples.json", "r") as file:
        file_info = json.load(file)
        public_files = []
        provided_files = []
        if include_publicly_available:
            public_files = [*map(lambda item: item["cleanFileName"], file_info["public"].values())]
        if include_provided:
            provided_files = [*map(lambda item: item["cleanFileName"], file_info["provided"].values())]
    return public_files + provided_files


def __poll_request(url: str) -> Response:
    wait_time: int = 1
    failed_request_count: int = 0
    while True:
        if failed_request_count >= NUM_FAILED_REQUESTS_ALLOWED:
            raise Exception("NCBI seems to be down at this time")
        response: Response = get(url, timeout=REQUEST_TIMEOUT)
        if response.status_code == 200:
            break
        if 400 <= response.status_code < 500 or response.status_code == 204:
            raise Exception("Something went wrong while fetching GEO data, please check the accessions")
        if response.status_code >= 500: 
            failed_request_count += 1
        print("Request failed...")
        time.sleep(wait_time)
        print("Retrying")
    return response

def __make_clean_file_name(title: str) -> str:
    for old, new in [(" ", ""), (":", "_"), ("+", "_"), (",", "_"), (";", "_")]:
        title = title.replace(old, new)
    return title

def get_macs_input():
    with open("config/samples.json") as file:
        macs_input = {}; samples = json.load(file)
        control_files: list[(str, str)] = []
        for _, sample in __flatten_dict(samples).items():
            name = f"{sample['mark']}_{sample['sample']}_rep{sample['replicate']}"
            if not sample['type'] == "input" and not name in macs_input:
                macs_input[name] = {"treatment": [], "control": []}
            if sample["type"] == "treatment":
                macs_input[name]["treatment"].append(sample["cleanFileName"])
                if not "peak_type" in macs_input[name]:
                    macs_input[name]["peak_type"] = sample["peak_type"]
                elif macs_input[name]["peak_type"] != sample["peak_type"]:
                    raise Exception("Mismatch between treatment sample peak types. Please check the sample sheet")
            elif sample["type"] == "control":
                control_files.append((name[2:], sample["cleanFileName"]))
    for name, file in control_files:
        for treatment_file in filter(lambda treatment_file: name in treatment_file, macs_input.keys()):
            macs_input[treatment_file]["control"].append(file)
    # Checks if a control sample has a treatment pair, if not ignores them
    macs_input = {key: value for key, value in macs_input.items() if value["treatment"] and value["control"]}
    return macs_input

def __flatten_dict(old_dict: dict) -> dict:
    new_dict = {}
    for key, value in old_dict.items():
        new_dict = new_dict | value
    return new_dict

if __name__ == "__main__":
    #make_sample_info()
    get_macs_input()