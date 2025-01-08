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

def join_read_files(runs: list, paired_end: bool):
    if paired_end:
        return [" ".join(list(map(lambda run: RESOURCES + f"/reads/{run}_1.fastq", runs))),
                " ".join(list(map(lambda run: RESOURCES + f"/reads/{run}_2.fastq", runs)))]
    return " ".join(list(map(lambda run: RESOURCES + f"/reads/{run}.fastq", runs))),

rule concatenate_runs_SE:
    input:
        lambda wildcards: expand(RESOURCES + "/reads/{run}.fastq",run=file_info["public"][wildcards.gsm]["runs"])
    output:
        RESOURCES + "/reads/{gsm}_{file_suffix}.fastq"
    wildcard_constraints:
        gsm = r"GSM[0-9]*",
    params:
        read_files = lambda wildcards: join_read_files(file_info["public"][wildcards.gsm]["runs"], False)
    log:
        LOGS + "/concatenate/{gsm}_{file_suffix}.log"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        cat {params.read_files} > {output} 
        """

rule concatenate_runs_PE:
    input:
        lambda wildcards: expand(RESOURCES + "/reads/{run}{read}.fastq", run=file_info["public"][wildcards.sample.split("_")[0]]["runs"], read=["_1", "_2"])
    output:
        read1 = RESOURCES + "/reads/{sample}_1.fastq",
        read2 = RESOURCES + "/reads/{sample}_2.fastq"
    params:
        read_files = lambda wildcards: join_read_files(file_info["public"][wildcards.sample.split("_")[0]]["runs"], True)
    log:
        LOGS + "/concatenate/{sample}.log"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        cat {params.read_files} > {output.read1} 
        cat {params.read_files} > {output.read2} 
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