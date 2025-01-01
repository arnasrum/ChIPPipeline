import sys
from sample_file_scripts import SampleFileScripts
sys.path.append("workflow/scripts")

file_info = SampleFileScripts.get_file_info()



read_ext = ["_1", "_2"] if str(config["paired_end"]).lower() == "true" else [""]
rule:
    name:
        f"download_SRR"
    output:
        temp(expand("resources/reads/{srr}{read}.fastq", read=read_ext, allow_missing=True))
    conda:
        "../envs/download.yml"
    wildcard_constraints:
        srr = r"SRR[0-9]*"
    log:
        "logs/fasterq-dump/{srr}"
    shell:
        '''
        exec > {log} 2>&1
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
            outputName = values["cleanFileName"],
            paired_end = config["paired_end"]
        log:
            f"logs/concatenate/{gsm}.log"
        shell:
            """
            exec > {log} 2>&1
            shopt -s nocasematch
            if [[ {params.paired_end} =~ true ]]; then
                cat {params.read1Files} > {params.outdir}/{params.outputName}_1.fastq
                cat {params.read2Files} > {params.outdir}/{params.outputName}_2.fastq
            else
                cat {params.readFiles} > {params.outdir}/{params.outputName}.fastq
            fi
            """


for key, value in file_info["provided"].items():
    reads = [file_info["provided"][key]["read1"]]
    if "read2" in file_info["provided"][key]: reads.append(file_info["provided"][key]["read2"])
    for read in reads:
        rule:
            name: f"link_{read['file_name']}"
            input:
                read["path"]
            output:
                f"resources/reads/{read['file_name']}{read['file_extension']}"
            params:
                pathToOriginal = read["path"],
                file_extension = read["file_extension"],
                file_name = read["file_name"]
            log:
                f"logs/link/{read['file_name']}.log"
            shell:
                '''
                exec > {log} 2>&1
                ln -s {input} {output} 
                '''

rule referenceGenome:
    output:
        "resources/genomes/{genome}.fa.gz"
    benchmark:
        "benchmarks/rsync/{genome}.benchmark.txt"
    conda:
        "../envs/download.yml"
    log:
        "logs/rsync/{genome}.log"
    shell:
        '''
        exec > {log} 2>&1
        rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz resources/genomes/
        '''