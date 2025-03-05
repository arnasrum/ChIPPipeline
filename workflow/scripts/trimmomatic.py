from snakemake.script import snakemake
from snakemake import shell
import sys

#if snakemake.threads:
    #command += f" -threads {int(snakemake.threads)}"

log = open(snakemake.log['log'], "w")
sys.stderr = log
sys.stdout = log

command = "trimmomatic"

def make_args():
    args = ""
    if snakemake.params.args:
        args += f" {snakemake.params['args']}"
        if snakemake.params['args'][-1] == " ":
            args += " "
    args += f" -threads {int(snakemake.threads)}"
    args += " -summary " + snakemake.log['summary']
    return args


if str(snakemake.params['paired_end']).lower() == 'true':
    command += " PE"
    command += make_args()
    command += f" {snakemake.input.samples[0]}"
    command += f" {snakemake.input.samples[1]}"
    command += f" {snakemake.output[1]}"
    command += f" {snakemake.output[0]}"
    command += f" {snakemake.output[3]}"
    command += f" {snakemake.output[2]}"
elif str(snakemake.params['paired_end']).lower() == 'false':
    command += " SE"
    command += make_args()
    command += f" {snakemake.input.samples[0]}"
    command += f" {snakemake.output[0]}"
else:
    raise Exception(f"Unknown parameter {snakemake.params['paired_end']}")

command += f" {snakemake.params['run_options']}"

log.close()
shell(f"exec > {snakemake.log['log']} 2>&1\n" + command)