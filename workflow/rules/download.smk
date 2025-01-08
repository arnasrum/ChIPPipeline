import sys
sys.path.append("workflow/scripts")
from sample_file_scripts import SampleFileScripts

sfs = SampleFileScripts(config)
file_info = SampleFileScripts.get_file_info(config["json_path"])

RESOURCES: str = config['resources_path']
LOGS: str = config['logs_path']
TEMP: str = config['temp_path']
BENCHMARKS: str = config['benchmarks_path']

rule fastq_dump_SE:
    output:
        temp(RESOURCES + "/reads/{srr}.fastq")
    params:
        path = f"{RESOURCES}/reads"
    conda:
        "../envs/download.yml"
    wildcard_constraints:
        srr = r"SRR[0-9]*"
    log:
        LOGS + "/fastq-dump/{srr}_SE.log"
    benchmark:
        BENCHMARKS + "/fastq-dump/{srr}_SE.log"
    resources:
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        fastq-dump -O {params.path} {wildcards.srr}
        '''

rule fastq_dump_PE:
    output:
        temp(RESOURCES + "/reads/{srr}_1.fastq"),
        temp(RESOURCES + "/reads/{srr}_2.fastq")
    params:
        path = f"{RESOURCES}/reads"
    conda:
        "../envs/download.yml"
    wildcard_constraints:
        srr = r"SRR[0-9]*"
    log:
        LOGS + "/fastq-dump/{srr}_PE.log"
    benchmark:
        BENCHMARKS + "/fastq-dump/{srr}_PE.log"
    resources:
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        fastq-dump -O {params.path} --split-3 {wildcards.srr}
        '''

read_ext = ["_1", "_2"] if sfs.is_paired_end() else [""]
for gsm, values in file_info["public"].items():
    rule:
        name:
            f"concatenate_{gsm}_SE"
        input:
            expand(RESOURCES + "/reads/{run}.fastq", run=values["runs"])
        output:
            RESOURCES + f"/reads/{values['file_name']}.fastq"
        params:
            path =  RESOURCES + "/reads",
            readFiles = " ".join(list(map(lambda run: RESOURCES + f"/reads/{run}.fastq", values["runs"]))),
            outputName = values["file_name"],
            paired_end = config["paired_end"]
        log:
            LOGS + f"/concatenate/{gsm}.log"
        resources:
            tmpdir=TEMP
        shell:
            """
            exec > {log} 2>&1
            cat {params.readFiles} > {params.path}/{params.outputName}.fastq
            """
    rule:
        name:
            f"concatenate_{gsm}_PE"
        input:
            expand(RESOURCES + "/reads/{run}{read}.fastq", run=values["runs"], read=["_1", "_2"])
        output:
            expand(RESOURCES + f"/reads/{values['file_name']}{{read}}.fastq", read=["_1", "_2"])
        params:
            path =  RESOURCES + "/reads",
            outputName= values["file_name"],
            paired_end= config["paired_end"],
            read1Files = " ".join(list(map(lambda run: RESOURCES + f"/reads/{run}_1.fastq", values["runs"]))),
            read2Files = " ".join(list(map(lambda run: RESOURCES + f"/reads/{run}_2.fastq", values["runs"])))
        log:
            LOGS + f"/concatenate/{gsm}.log"
        resources:
            tmpdir=TEMP
        shell:
            """
            exec > {log} 2>&1
            cat {params.read1Files} > {params.path}/{params.outputName}_1.fastq
            cat {params.read2Files} > {params.path}/{params.outputName}_2.fastq
            """

for key, value in file_info["provided"].items():
    rule:
        name: f"link_{value['file_name']}_SE"
        input:
            [file_info["provided"][key]["read1"]["path"].rstrip(".gz")]
        output:
            RESOURCES + f"/reads/{value['file_name']}.fastq"
        log:
            LOGS + f"/link/{value['file_name']}.log"
        resources:
            tmpdir=TEMP
        shell:
            '''
            exec > {log} 2>&1
            ln -s {input} {output} 
            '''
    rule:
        name: f"symlink_{value['file_name']}_PE"
        input:
            files = [file_info["provided"][key]["read1"]["path"].rstrip(".gz"), file_info["provided"][key]["read2"]["path"].rstrip(".gz")]
        output:
            out_files = [RESOURCES + f"/reads/{value['file_name']}_1.fastq", RESOURCES + f"/reads/{value['file_name']}_2.fastq"]
        log:
             LOGS + f"/link/{value['file_name']}.log"
        resources:
            tmpdir=TEMP
        shell:
            '''
            exec > {log} 2>&1
            ln -sr {input.files[0]} {output.out_files[0]} 
            ln -sr {input.files[1]} {output.out_files[1]} 
            '''

rule referenceGenome:
    output:
        RESOURCES + "/genomes/{genome}.fa.gz"
    conda:
        "../envs/download.yml"
    log:
        LOGS + "/rsync/{genome}.log"
    benchmark:
        BENCHMARKS + "/rsync/{genome}.benchmark.txt"
    resources:
        tmpdir=TEMP
    shell:
        '''
        exec > {log} 2>&1
        rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz {output}
        '''