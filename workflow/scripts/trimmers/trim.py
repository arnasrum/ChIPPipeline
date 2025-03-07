from snakemake.script import snakemake
from snakemake import shell
from workflow.scripts.trimmers.fastp import Fastp
from workflow.scripts.trimmers.trimmomatic import Trimmomatic
from workflow.scripts.trimmers.cutadapt import Cutadapt
from workflow.scripts.trimmers.trim_galore import TrimGalore
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
if str(snakemake.config['trimmer']).lower() == "fastp":
    trimmer = Fastp()
elif str(snakemake.config['trimmer']).lower() == "trimmomatic":
    trimmer = Trimmomatic()
    trimmer.run_options = snakemake.config["trimmomatic"]["run_options"]
elif str(snakemake.config['trimmer']).lower() == "cutadapt":
    trimmer = Cutadapt()
elif str(snakemake.config['trimmer']).lower() == "trim_galore":
    trimmer = TrimGalore()
else:
    raise Exception(f"Trimmer {snakemake.config['trimmer']} not supported")

read1 = snakemake.input[0]
read2 = snakemake.input[1] if snakemake.input[1] is not None else None
threads = int(snakemake.threads)
args = snakemake.config[snakemake.config['trimmer']]['args'] if snakemake.params.args is not None else None
output_path = str(os.path.commonpath(snakemake.output))
command = trimmer.trim(output_path, read1, read2, threads, args)
shell(f"({command}) {log}")