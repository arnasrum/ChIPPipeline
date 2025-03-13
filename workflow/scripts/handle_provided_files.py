from snakemake.script import snakemake
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command = ""
if snakemake.input[0][-3:] == ".gz":
    command = f"ln -sr {snakemake.input[0]} {snakemake.output[0]}"
else:
    command = f"gzip -kfc9 {snakemake.input[0]} > {snakemake.output[0]}"
if len(snakemake.output) == 2:
    if snakemake.input[1][-3:] == ".gz":
        command += f"\nln -sr {snakemake.input[1]} {snakemake.output[1]}"
    else:
        command += f"\ngzip -kfc9 {snakemake.input[1]} > {snakemake.output[1]}"

shell(f"({command}) {log}")