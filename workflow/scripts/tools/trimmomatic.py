from snakemake.script import snakemake
from snakemake import shell
import os


log = snakemake.log_fmt_shell(stdout=True, stderr=True)

read2 = None
file_name = snakemake.wildcards['sample']
output_path = os.path.dirname(snakemake.output[0])
read_extensions = [".fastq"]
if len(snakemake.input) == 2:
    read_extensions = ["_1.fastq", "_2.fastq"]
    read2 = snakemake.input[1]


command = "trimmomatic"
if read2:
    command += " PE"
else:
    command += " SE"
command += f" -threads {snakemake.threads}"
command += f" -summary {output_path}/{file_name}_summary.txt"
if snakemake.params.args:
    command += f" {snakemake.params.args}"
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
command += f" {snakemake.params.run_options}"
if read2:
    command += f"\nrm {unpaired1}"
    command += f"\nrm {unpaired2}"

shell(f"({command}) {log}")