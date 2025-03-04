from snakemake.script import snakemake
from snakemake import shell
import sys

#if snakemake.threads:
    #command += f" -threads {int(snakemake.threads)}"

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

command = "trimmomatic"
if snakemake.params.args:
    command += f" {snakemake.params['args']}"
    if snakemake.params['args'][-1] == " ":
        command += " "
if str(snakemake.params['paired_end']).lower() == 'true':
    command += " PE"
    command += f" -threads {int(snakemake.threads)}"
    command += f" {snakemake.input.samples[0]}"
    command += f" {snakemake.input.samples[1]}"
    command += f" {snakemake.output[0]}"
    command += f" {snakemake.output[1]}"
elif str(snakemake.params['paired_end']).lower() == 'false':
    command += " SE"
    command += f" -threads {int(snakemake.threads)}"
    command += f" {snakemake.input.samples[0]}"
    command += f" {snakemake.output[0]}"
command += f" {snakemake.params['run_options']}"

log.close()
shell(f"exec > {snakemake.log[0]} 2>&1\n" + command)