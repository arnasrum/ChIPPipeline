from snakemake.script import snakemake
from snakemake import shell
import os
from bowtie2 import Bowtie2
from bwa_mem2 import BwaMem2
from workflow.scripts.aligners.star import STAR

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.config['aligner'] == "bowtie2":
    aligner = Bowtie2()
elif snakemake.config['aligner'] == "bwa-mem2":
    aligner = BwaMem2()
elif str(snakemake.config['aligner']).lower() == "star":
    aligner = STAR()
else:
   raise Exception(f"Unknown aligner {snakemake.config['aligner']}")
prefix = str(os.path.commonprefix(snakemake.output).rstrip("."))
shell(f"({aligner.build_index(snakemake.input, prefix, snakemake.threads, snakemake.params.args)}) {log}")
