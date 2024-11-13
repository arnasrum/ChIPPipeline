import sys
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
        conda:
            "../envs/download.yml"
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


reads = ["1", "2"] if config["paired_end"] else ["1"]
for key, value in fileInfo["provided"].items():
    rule:
        name: f"link_{value["cleanFileName"]}"
        input:
            expand("{path}{fileName}_{num}{ext}", fileName=value["cleanFileName"], path=value["path"], num=reads, ext=value["fileExtension"])
        output:
            expand("resources/reads/{fileName}_{num}{ext}", fileName=value["cleanFileName"], num=reads, ext=value["fileExtension"])
        params:
            paired_end = config["paired_end"],
            pathToOriginal = f"{value['path']}{value['cleanFileName']}",
            fileExt = value["fileExtension"],
            cleanFileName = value["cleanFileName"],
        shell:
            '''
                if [[ {params.paired_end} ]]; then
                    ln {params.pathToOriginal}_1{params.fileExt} resources/reads/{params.cleanFileName}_1{params.fileExt}
                    ln {params.pathToOriginal}_2{params.fileExt} resources/reads/{params.cleanFileName}_2{params.fileExt}
                elif [[ {params.paired_end} ]]; then
                    ln {params.pathToOriginal}_1{params.fileExt} resources/reads/{params.cleanFileName}_1{params.fileExt}
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