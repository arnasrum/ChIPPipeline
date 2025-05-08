from snakemake.script import snakemake
from snakemake import shell
from os import path
from uuid import uuid4

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command = "STAR"

command += f" --outTmpDir {snakemake.resources['tmpdir']}/{str(uuid4())}"
output_path = path.dirname(snakemake.output[0]) + f"/{snakemake.wildcards['sample']}_"
command += f" --outFileNamePrefix {output_path}"
command += f" --readFilesType Fastx"
command += f" --runThreadN {snakemake.threads}"
genome_dir = str(path.commonpath(snakemake.input['genome_index']))
command += f" --genomeDir {genome_dir}"
command += f" --readFilesIn {snakemake.input['reads']}"
if snakemake.params['args']:
    command += f" {snakemake.params['args']}"
command += f" --outSAMtype BAM Unsorted"
command += f" --outSAMattrRGline ID:{snakemake.wildcards['sample']} SM:{snakemake.wildcards['sample']}"
command += f"\nmv {snakemake.output[0].replace('.bam', '')}_Aligned.out.bam {snakemake.output[0]}"

shell(f"({command}) {log}")
