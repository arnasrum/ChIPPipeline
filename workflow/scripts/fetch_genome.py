from snakemake.script import snakemake
from snakemake import shell
import pandas as pd
from  os import path

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command =  ""
with open(snakemake.config['sample_sheet']) as sample_sheet:
    samples = pd.read_csv(sample_sheet, header=0)['genome'].tolist()
file = next(filter(lambda genome: snakemake.params['genome'] in genome and path.isfile(genome), samples), None)
if file:
    if file[-3:] == ".gz":
        command += "echo 'gzipped genome detected; symlinking'"
        command += f"\nln -sr {file} {snakemake.output[0]}"
    else :
        command += "echo 'non-gzipped genome detected; gzipping'"
        command += f"\ngzip -kfc {file} > {snakemake.output[0]}"
else:
    command += "echo 'Defined genome is not a file; trying to download'"
    command += f"\nrsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/{snakemake.params['genome']}/bigZips/{snakemake.params['genome']}.fa.gz {snakemake.output[0]}"
shell(f"({command}) {log}")
