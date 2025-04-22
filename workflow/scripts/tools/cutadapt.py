from snakemake.script import snakemake
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
command = "cutadapt"
command += f" -j {snakemake.threads}"
command += f" -o {snakemake.output[0]}"
if len(snakemake.output) == 2:
    command += f" -p {snakemake.output[1]}"
if snakemake.config['cutadapt']['args']:
    command += f" {snakemake.config['cutadapt']['args']}"
command += f" {snakemake.input[0]}"
if len(snakemake.input) == 2:
    command += f" {snakemake.input[1]}"

shell(f"({command}) {log}")