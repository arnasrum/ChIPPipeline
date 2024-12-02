import sys
from sample_file_scripts import make_sample_info
sys.path.append("workflow/scripts")

file_info = make_sample_info()



read_ext = ["_1", "_2"] if str(config["paired_end"]).lower() == "true" else [""]
rule:
    name:
        f"download_SRR"
    output:
        temp(expand("resources/reads/{srr}{read}.fastq", read=read_ext, allow_missing=True))
    conda:
        "../envs/download.yml"
    shell:
        '''
        fasterq-dump --temp temp -O resources/reads -p {wildcards.srr}
        '''


for gsm, values in file_info["public"].items():
    rule:
        name:
            f"concatenate_{gsm}"
        input:
            expand("resources/reads/{run}{readNum}.fastq", run=values["runs"], readNum=read_ext)
        output:
            expand(f"resources/reads/{values['cleanFileName']}{{read}}.fastq", read=read_ext)
        params:
            outdir = "resources/reads",
            readFiles = " ".join(list(map(lambda run: f"resources/reads/{run}.fastq", values["runs"]))),
            read1Files = " ".join(list(map(lambda run: f"resources/reads/{run}_1.fastq", values["runs"]))),
            read2Files = " ".join(list(map(lambda run: f"resources/reads/{run}_2.fastq", values["runs"]))),
            outputName = values["cleanFileName"] ,
            paired_end = config["paired_end"]
        shell:
            """
            shopt -s nocasematch
            if [[ {params.paired_end} =~ true ]]; then
                cat {params.read1Files} > {params.outdir}/{params.outputName}_1.fastq
                cat {params.read2Files} > {params.outdir}/{params.outputName}_2.fastq
            else
                cat {params.readFiles} > {params.outdir}/{params.outputName}.fastq
            fi
            """


reads = ["_1", "_2"] if config["paired_end"] else [""]
for key, value in file_info["provided"].items():
    rule:
        name: f"link_{value['cleanFileName']}"
        input:
            expand("{path}{fileName}_{num}{ext}", fileName=value["cleanFileName"], path=value["path"], num=reads, ext=value["fileExtension"])
        output:
            expand("resources/reads/{fileName}{num}{ext}", fileName=value["cleanFileName"], num=reads, ext=value["fileExtension"])
        params:
            paired_end = config["paired_end"],
            pathToOriginal = f"{value['path']}{value['cleanFileName']}",
            fileExt = value["fileExtension"],
            cleanFileName = value["cleanFileName"]
        shell:
            '''
                shopt -s nocasematch
                if [[ {params.paired_end} =~ true ]]; then
                    ln {params.pathToOriginal}_1{params.fileExt} resources/reads/{params.cleanFileName}_1{params.fileExt}
                    ln {params.pathToOriginal}_2{params.fileExt} resources/reads/{params.cleanFileName}_2{params.fileExt}
                else
                    ln {params.pathToOriginal}_1{params.fileExt} resources/reads/{params.cleanFileName}{params.fileExt}
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