from snakemake.script import snakemake
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
command = "fastp"
command += f" -w {snakemake.threads}"
command += f" -i {snakemake.input[0]}"
command += f" -o {snakemake.output[0]}"
if len(snakemake.input) == 2:
    command += f" -I {snakemake.input[1]}"
    command += f" -O {snakemake.output[1]}"
command += f" -j {snakemake.output[0].split('.')[0]}.json"
command += f" -h {snakemake.output[0].split('.')[0]}.html"
if snakemake.config["fastp"]["args"]:
    command += " " + snakemake.config["fastp"]["args"]
shell(f"({command}) {log}")