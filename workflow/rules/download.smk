import sys
import pandas as pd
import numpy as np
from sampleFileScripts import makeSampleInfo 
sys.path.append("workflow/scripts")

fileInfo = makeSampleInfo()

for srr in [run for value in fileInfo["public"].values() for run in value["runs"]]:
    rule:
        name:
            f"download_{srr}"
        output:
            temp(f"resources/reads/{srr}_1.fastq"),
            temp(f"resources/reads/{srr}_2.fastq")
        params:
            srr = srr
        shell:
            '''
            fasterq-dump --temp temp -O resources/reads -p {params.srr}
            '''


for gsm, values in fileInfo["public"].items():
    rule:
        name:
            f"download_{gsm}"
        input:
            expand("resources/reads/{run}_{readNum}.fastq", run=values["runs"], readNum=[1, 2])
        output:
            f"resources/reads/{values["cleanFileName"]}_1.fastq",
            f"resources/reads/{values["cleanFileName"]}_2.fastq"
        params:
            outdir = "resources/reads",
            read1Files = " ".join(list(map(lambda run: f"resources/reads/{run}_1.fastq", values["runs"]))),
            read2Files = " ".join(list(map(lambda run: f"resources/reads/{run}_2.fastq", values["runs"]))),
            outputName = values["cleanFileName"] 
        shell:
            """
            cat {params.read1Files} > {params.outdir}/{params.outputName}_1.fastq
            cat {params.read2Files} > {params.outdir}/{params.outputName}_2.fastq
            """


reads = ["1", "2"] if config["libraryStrategy"] == "paired" else ["1"]
for key, value in fileInfo["provided"].items():
    rule:
        name: f"link_{value["cleanFileName"]}"
        input:
            expand("{path}_{num}.{ext}", path=value["path"], num=reads, ext=value["fileExtension"]) 
        output:
            expand("resources/reads/{fileName}_{num}.fastq", fileName=value["cleanFileName"], num=reads) 
        params:
            libraryStrategy = config["libraryStrategy"],
            pathToOriginal = value["path"],
            fileExt = value["fileExtension"],
            cleanFileName = value["cleanFileName"]
        shell:
            '''
                if [[ {params.libraryStrategy} == "paired" ]]; then
                    ln {params.pathToOriginal}_1.{params.fileExt} resources/reads/{params.cleanFileName}_1.fastq
                    ln {params.pathToOriginal}_2.{params.fileExt} resources/reads/{params.cleanFileName}_2.fastq
                elif [[ {params.libraryStrategy} == "single" ]]; then
                    ln {params.pathToOriginal}_1.{params.fileExt} resources/reads/{params.cleanFileName}_1.fastq
                fi
            '''

rule referenceGenome:
    output:
        "resources/genomes/{genome}.fa.gz"
    benchmark:
        "benchmarks/rsync/{genome}.benchmark.txt"
    shell:
        '''
        rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz resources/genomes/
        '''