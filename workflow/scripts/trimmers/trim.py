from snakemake.script import snakemake
from snakemake import shell
from fastp import Fastp
from trimmomatic import Trimmomatic
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if str(snakemake.config['trimmer']).lower() == "fastp":
    trimmer = Fastp()
elif str(snakemake.config['trimmer']).lower() == "trimmomatic":
    trimmer = Trimmomatic()
    trimmer.run_options = snakemake.config["trimmomatic"]["run_options"]
else:
    raise Exception(f"Trimmer {snakemake.config['trimmer']} not supported")

read1 = snakemake.input[0]
read2 = snakemake.input[1] if snakemake.input[1] is not None else None
threads = int(snakemake.threads)
args = snakemake.params.args if snakemake.params.args is not None else None
output_path = str(os.path.commonpath(snakemake.output))
command = trimmer.trim(output_path, read1, read2, threads, args)
shell(f"({command}) {log}")