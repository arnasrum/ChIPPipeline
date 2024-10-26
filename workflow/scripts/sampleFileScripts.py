from requests import get, Response
from xml.etree import ElementTree
from typing import Any
import pandas as pd
import pathlib
import json
import os
import re

REQUEST_TIMEOUT = 30

def makeSampleInfo(sampleSheet:str="config/samples.csv") -> dict[str:dict]:

    #if os.path.isfile("config/samples.json"):
        #with open("config/samples.json", "r") as file:
            #return json.load(file)

    pattern = re.compile(r"^GSM[0-9]*$")
    inputSamples = set()
    with open(sampleSheet, "r") as file:
        for index, row in pd.read_csv(file).iterrows():
            for gsm in row.values[2:]:
                inputSamples.add(gsm)
    sampleInfo = {}
    sampleInfo["public"] = {}; sampleInfo["provided"] = {}
    providedInfo: dict[str:Any] = sampleInfo["provided"]
    gsmAccessions = set()

    gsmAccessions = set(filter(lambda item: pattern.match(item), inputSamples))
    # Handle provided files
    # gsm should be renamed
    for gsm in filter(lambda item: not pattern.match(item), inputSamples):
        path = pathlib.Path(gsm)
        fileExtention = "".join(path.suffixes)
        fileName = path.name.split(fileExtention)[0]
        providedInfo[fileName] = {}
        providedInfo[fileName]["cleanFileName"] = fileName
        providedInfo[fileName]["fileExtension"] = fileExtention 
        providedInfo[fileName]["path"] = gsm.split(fileName + fileExtention)[0]
    # Handle publicly available files
    sampleInfo["public"] = {key:value for key, value in getMetaData(getSraAccessions(gsmAccessions).values()).items()}

    with open(f"config/samples.json", "w") as outfile:
        outfile.write(json.dumps(sampleInfo, indent=4))
    return sampleInfo

def getAllSampleFilePaths(includeDirectories=True) -> list[str]:
    directory = "resources/reads/" if includeDirectories else ""
    filePaths = []
    with open("config/samples.json", "r") as file:
        data: dict[str:dict] = json.load(file)
        for type in data:
            for fileInfo in data[type].values():
                path = f"{directory}{fileInfo["cleanFileName"]}"
                filePaths.append(f"{path}_1.fastq")
                filePaths.append(f"{path}_2.fastq")

    return filePaths


def getSraAccessions(geoAccessions: list[str]) -> dict[str:str]:
    '''
        Given a list of GEO accessions retrieve their related SRA accessions 
    '''
    if len(geoAccessions) == 0: return dict() 
    enterezUrl: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term="
    enterezUrl += "+OR+".join(geoAccessions)
    print(enterezUrl)
    ### Handle 5** and 4** status codes
    response: Response = get(enterezUrl, timeout=REQUEST_TIMEOUT)
    xmlResponse: ElementTree.Element = ElementTree.fromstring(response.content)
    idList = list(filter(lambda item: int(item) > 299999999, [id.text for id in xmlResponse[3]]))
    enterezUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&id=" + ",".join(idList)
    print(enterezUrl)
    ### Handle 5** and 4** status codes
    response = get(enterezUrl, timeout=REQUEST_TIMEOUT)
    responseText: str = response.content.decode().lstrip("\n").rstrip("\n")
    gsmToSraMap: dict[str:str] = {} 
    responseText = responseText.replace("\n\n", "\t;:.,").replace("\n", "")
    for match in responseText.split("\t;:.,"):
        subResponse: str = match 
        accessionPosition: re.Match = re.search(r"Accession:\sGSM[0-9]*", subResponse)
        sraAccessionPosition: re.Match = re.search(r"SRX[0-9]*", subResponse)
        accession: str = subResponse[accessionPosition.span()[0]: accessionPosition.span()[1]]
        sraAccession: str = subResponse[sraAccessionPosition.span()[0]: sraAccessionPosition.span()[1]]
        gsmToSraMap[accession] = sraAccession
    return gsmToSraMap

def getMetaData(sraAccessions: list[str]) -> dict[str: dict]: 
    '''
        Fetch metadata for SRA samples 
    '''
    if len(sraAccessions) == 0:
        return dict()
    runAccessions: list[str] = []
    enterezUrl: str = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={",".join(sraAccessions)}"
    print(enterezUrl)
    ### Handle 5** and 4** status codes
    response: Response = get(enterezUrl, timeout=REQUEST_TIMEOUT)
    root: ElementTree.Element = ElementTree.fromstring(response.content)
    metaData: dict[str: dict] = {}
    for node in root:
        runAccessions: list[str] = []
        geoAccession = node[4].attrib["alias"]
        metaData[geoAccession] = {"cleanFileName": __makeCleanFileName(node[0][1].text)}
        for run in node[6]:
            runAccessions.append(run.attrib["accession"])
        metaData[geoAccession]["runs"] = runAccessions
    #response = xmltodict.parse(result.content)
    return metaData

def getFileNames(includeProvided: bool = True, includePubliclyAvailiable: bool = True) -> list[str]:
    """
        Returns a list of file names 
    """
    if not includeProvided and not includePubliclyAvailiable:
        return list()

    if not os.path.isfile("config/samples.json"):
        makeSampleInfo()
    #makeSampleInfo()
    with open("config/samples.json", "r") as file:
        fileInfo = json.load(file)
        publicFiles = []
        providedFiles = []
        if includePubliclyAvailiable:
            publicFiles = [*map(lambda item: item["cleanFileName"], fileInfo["public"].values())]
        if includeProvided:
            providedFiles = [*map(lambda item: item["cleanFileName"], fileInfo["provided"].values())]
    return publicFiles + providedFiles


def __makeCleanFileName(title: str) -> str:
    for old, new in [(" ", ""), (":", "_"), ("+", "_"), (",", "_"), (";", "_")]:
        title = title.replace(old, new)
    return title

if __name__ == "__main__":
    makeSampleInfo()
    #print(getAllFileNames())
    #print(getFileNames(includeProvided=True, includePubliclyAvailiable=True))