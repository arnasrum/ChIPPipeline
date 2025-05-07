from snakemake.script import snakemake
from snakemake import shell
import os


log = snakemake.log_fmt_shell(stdout=True, stderr=True)

read2 = None
file_name = snakemake.wildcards['sample']
output_path = os.path.dirname(snakemake.output[0])
read_extensions = [".fastq.gz"]
if len(snakemake.input) == 2:
    read_extensions = ["_1.fastq.gz", "_2.fastq.gz"]
    read2 = snakemake.input[1]

command = "trimmomatic"
if read2:
    command += " PE"
else:
    command += " SE"
command += f" -threads {snakemake.threads}"
command += f" -summary {output_path}/{file_name}_summary.txt"
if snakemake.config['trimmomatic']['args']:
    command += f" {snakemake.config['trimmomatic']['args']}"
command += f" {snakemake.input[0]}"
if read2:
    command += f" {read2}"
command += f" {output_path}/{file_name}{read_extensions[0]}"
if read2:
    unpaired1 = f"{output_path}/{file_name}_unpaired{read_extensions[0]}"
    command += f" {unpaired1}"
    command += f" {output_path}/{file_name}{read_extensions[1]}"
    unpaired2 = f"{output_path}/{file_name}_unpaired{read_extensions[1]}"
    command += f" {unpaired2}"
if snakemake.config['trimmomatic']['trimming_steps']:
    command += f" {snakemake.config['trimmomatic']['run_options']}"
elif read2:
    command += f" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36"
else:
    command += f" ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
if read2:
    command += f"&& rm {unpaired1}"
    command += f"&& rm {unpaired2}"

shell(f"({command}) {log}")