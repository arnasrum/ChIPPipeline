from snakemake.script import snakemake
from snakemake import shell
from os import path

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command =  ""
path: str | None = snakemake.params['path']
if path:
    if path.endswith(".gz"):
        command += "echo 'gzipped genome detected; symlinking'"
        command += f"\nln -sr {path} {snakemake.output[0]}"
    else :
        command += "echo 'non-gzipped genome detected; gzipping'"
        command += f"\ngzip -kfc {path} > {snakemake.output[0]}"
else:
    command += "echo 'Defined genome is not a file; trying to download'"
    command += f"\nrsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{snakemake.wildcards['genome']}/bigZips/{snakemake.wildcards['genome']}.fa.gz {snakemake.output[0]}"
shell(f"({command}) {log}")
