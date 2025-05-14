from snakemake.script import snakemake
from snakemake import shell
from os import path

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command = "bwa mem"
command += f" -t {snakemake.threads}"
if snakemake.params['args']:
    command += f" {snakemake.params['args']}"
command += f" -R '@RG\\tID:{snakemake.wildcards.sample}\\tSM:{snakemake.wildcards.sample}'"
command += f" {path.commonprefix(snakemake.input['genome_index'])}".rstrip(".")
command += " " + " ".join(snakemake.input['reads'])
command += f" | samtools view -o {snakemake.output[0]} - "
shell(f"(echo 'Executing: {command}'\n{command}) {log}")