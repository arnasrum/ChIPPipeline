from snakemake.script import snakemake
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command = ""
if snakemake.input[0][-3:] == ".gz":
    command = f"gzip -dkf {snakemake.input[0]}"
    command += f"\nmv {snakemake.input[0].replace(".gz", "")} {snakemake.output[0]}"
else:
    command = f"ln -sr {snakemake.input[0]} {snakemake.output[0]}"
if str(snakemake.config['paired_end']).lower() == "true":
    if snakemake.input[1][-3:] == ".gz":
        command += f"\ngzip -dkf {snakemake.input[1]}"
        command += f"\nmv {snakemake.input[1].replace(".gz", "")} {snakemake.output[1]}"
    else:
        command += f"\nln -sr {snakemake.input[1]} {snakemake.output[1]}"

shell(f"(echo 'Running: {command}'\n{command}) {log}")