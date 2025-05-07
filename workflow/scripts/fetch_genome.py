from snakemake.script import snakemake
from snakemake import shell
from os import path

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command =  ""
file = next(filter(lambda genome: snakemake.wildcards['genome'] in genome and path.isfile(genome), snakemake.params['samples']), None)
if file:
    if file[-3:] == ".gz":
        command += "echo 'gzipped genome detected; symlinking'"
        command += f"\nln -sr {file} {snakemake.output[0]}"
    else :
        command += "echo 'non-gzipped genome detected; gzipping'"
        command += f"\ngzip -kfc {file} > {snakemake.output[0]}"
else:
    command += "echo 'Defined genome is not a file; trying to download'"
    command += f"\nrsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{snakemake.wildcards['genome']}/bigZips/{snakemake.wildcards['genome']}.fa.gz {snakemake.output[0]}"
shell(f"({command}) {log}")
