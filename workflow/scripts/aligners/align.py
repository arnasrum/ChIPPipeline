from snakemake.script import snakemake
from snakemake import shell
from workflow.scripts.aligners.bowtie2 import Bowtie2
from workflow.scripts.aligners.bwa_mem2 import BwaMem2
from workflow.scripts.aligners.star import STAR
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.config['aligner'] == "bowtie2":
    aligner = Bowtie2()
elif snakemake.config['aligner'] == "bwa-mem2":
    aligner = BwaMem2()
elif str(snakemake.config['aligner']).lower() == "star":
    aligner = STAR()
    aligner.output_prefix = os.path.dirname(snakemake.output[0])
else:
    raise Exception(f"Unknown aligner {snakemake.config['aligner']}")

index = str(os.path.commonprefix(snakemake.input['index']).rstrip("."))
read1 = snakemake.input['reads'][0]
read2 = snakemake.input['reads'][1] if len(snakemake.input['reads']) == 2 else None
command = aligner.align(index, read1, read2, int(snakemake.threads), str(snakemake.params.args))
command += f" | samtools sort -@ {snakemake.threads} -o {snakemake.output} -"
shell(f"({command}) {log}")
