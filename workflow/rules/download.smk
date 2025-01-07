import sys
from sample_file_scripts import SampleFileScripts
sys.path.append("workflow/scripts")

sfs = SampleFileScripts(config)
file_info = SampleFileScripts.get_file_info()

rule fastq_dump_SE:
    output:
        temp("resources/reads/{srr}.fastq")
    conda:
        "../envs/download.yml"
    wildcard_constraints:
        srr = r"SRR[0-9]*"
    log:
        "logs/fastq-dump/{srr}_SE.log"
    shell:
        '''
        exec > {log} 2>&1
        fastq-dump -O resources/reads {wildcards.srr}
        '''

rule fastq_dump_PE:
    output:
        temp("resources/reads/{srr}_1.fastq"),
        temp("resources/reads/{srr}_2.fastq")
    conda:
        "../envs/download.yml"
    wildcard_constraints:
        srr = r"SRR[0-9]*"
    log:
        "logs/fastq-dump/{srr}_PE.log"
    shell:
        '''
        exec > {log} 2>&1
        fastq-dump -O resources/reads --split-3 {wildcards.srr}
        '''


for gsm, values in file_info["public"].items():
    rule:
        name:
            f"concatenate_{gsm}"
        input:
            expand("resources/reads/{run}{readNum}.fastq", run=values["runs"], readNum=read_ext)
        output:
            expand(f"resources/reads/{values['file_name']}{{read}}.fastq", read=read_ext)
        params:
            outdir = "resources/reads",
            readFiles = " ".join(list(map(lambda run: f"resources/reads/{run}.fastq", values["runs"]))),
            read1Files = " ".join(list(map(lambda run: f"resources/reads/{run}_1.fastq", values["runs"]))),
            read2Files = " ".join(list(map(lambda run: f"resources/reads/{run}_2.fastq", values["runs"]))),
            outputName = values["file_name"],
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
    rule:
        name: f"link_{value['file_name']}_SE"
        input:
            [file_info["provided"][key]["read1"]["path"].rstrip(".gz")]
        output:
            [f"resources/reads/{value['file_name']}.fastq"]
        log:
            f"logs/link/{value['file_name']}.log"
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
            out_files = [f"resources/reads/{value['file_name']}_1.fastq",f"resources/reads/{value['file_name']}_2.fastq"]
        log:
            f"logs/link/{value['file_name']}.log"
        shell:
            '''
            exec > {log} 2>&1
            ln -sr {input.files[0]} {output.out_files[0]} 
            ln -sr {input.files[1]} {output.out_files[1]} 
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