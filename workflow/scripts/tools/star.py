from snakemake.script import snakemake
from snakemake import shell
from os import path
from uuid import uuid4

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command = "STAR"

uuid = str(uuid4())
command += f" --outTmpDir {snakemake.resources['tmpdir']}/{uuid}"
output_path = path.dirname(snakemake.output[0]) + snakemake.wildcards['sample'] + "_"
command += f" --outFileNamePrefix {output_path}"
command += f" --readFilesType Fastx"
command += f" --runThreadN {snakemake.threads}"
genome_dir = str(path.commonpath(snakemake.input['genome_index']))
command += f" --genomeDir {genome_dir}"
command += f" --readFilesIn {snakemake.input['reads']}"
if snakemake.params['args']:
    command += f" {snakemake.params['args']}"
command += f" --outStd BAM"
command += f" > {snakemake.output[0]}"
shell(f"({command}) {log}")



"""
exec > {log} 2>&1
uuid=$(python3 -c "import uuid; print(uuid.uuid4())")
STAR --outTmpDir "{resources.tmpdir}/STAR-${{uuid}}" --outFileNamePrefix {params.output_path}/{wildcards.sample}_ --readFilesType Fastx --runThreadN {threads} --genomeDir {params.index_path} --readFilesIn {input.reads} {params.args} --outStd SAM | samtools sort -@ {threads} -o {output} -
"""

