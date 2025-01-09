from requests import get, Response
from xml.etree import ElementTree
import time
import json
import re
import os

REQUEST_TIMEOUT: int = 30
NUM_FAILED_REQUESTS_ALLOWED: int = 5


def get_sra_accessions(geo_accessions: set[str]) -> dict[str:str]:
    """
        Given a list of GEO accessions retrieve their related SRA accessions
    """
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
    """
        Fetch metadata for SRA samples
    """
    if len(sra_accessions) == 0:
        return dict()
    enterez_url: str = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={','.join(sra_accessions)}"
    print(enterez_url)
    response = __poll_request(enterez_url)
    root: ElementTree.Element = ElementTree.fromstring(response.content)
    meta_data: dict[str: dict] = {}
    for node in root:
        run_accessions: list[str] = []
        geo_accession = node[4].attrib["alias"]
        meta_data[geo_accession] = {"file_name": __make_clean_file_name(node[0][1].text)}
        for run in node[6]:
            run_accessions.append(run.attrib["accession"])
        meta_data[geo_accession]["runs"] = run_accessions
    #response = xmltodict.parse(result.content)
    return meta_data




def __poll_request(url: str) -> Response:
    wait_time: int = 1
    failed_request_count: int = 0
    while True:
        if failed_request_count >= NUM_FAILED_REQUESTS_ALLOWED:
            raise Exception("NCBI seems to be down at this time")
        response: Response = get(url, timeout=REQUEST_TIMEOUT)
        if response.status_code == 200 and len(response.content) > 299:
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
    for old, new in [(" ", ""), (":", "_"), ("+", "_"), (",", "_"), (";", "_"), (".", "_")]:
        title = title.replace(old, new)
    return title.rstrip("_1").rstrip("_2")
