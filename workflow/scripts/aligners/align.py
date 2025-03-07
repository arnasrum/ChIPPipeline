from snakemake.script import snakemake
from snakemake import shell
from bowtie2 import Bowtie2
from bwa_mem2 import BwaMem2
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.config['aligner'] == "bowtie2":
    aligner = Bowtie2()
elif snakemake.config['aligner'] == "bwa-mem2":
    aligner = BwaMem2()
else:
    raise Exception(f"Unknown aligner {snakemake.config['aligner']}")

prefix = str(os.path.commonprefix(snakemake.input['index']).rstrip("."))
read1 = snakemake.input['reads'][0]
read2 = snakemake.input['reads'][1] if len(snakemake.input['reads']) == 2 else None
shell(f"({aligner.align(prefix, read1, read2, snakemake.threads, snakemake.params.args)}) {log}")
