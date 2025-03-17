from snakemake.script import snakemake
from snakemake import shell
from os import path


log = snakemake.log_fmt_shell(stdout=True, stderr=True)
command = "bowtie2 -q"
command += f" --rg-id {snakemake.wildcards['sample']}"
command += f" --rg SM:{snakemake.wildcards['sample']}"

command += f" --threads {snakemake.threads}"
if snakemake.params['args']:
    command += f" {snakemake.params['args']}"

command += f" -x {path.commonprefix(snakemake.input['genome_index']).rstrip('.')}"
if len(snakemake.input['reads']) == 2:
    command += f" -1 {snakemake.input['reads'][0]} -2 {snakemake.input['reads'][1]}"
elif len(snakemake.input['reads']) == 1:
    command += f" -U {snakemake.input['reads'][0]}"

shell(f"({command} | samtools sort -@ {snakemake.threads} -o {snakemake.output[0]} -) {log}")