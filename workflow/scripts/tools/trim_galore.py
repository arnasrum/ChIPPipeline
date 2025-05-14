from snakemake.script import snakemake
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
command = "trim_galore"
command += f" --basename {snakemake.wildcards.sample}"
command += f" -j {snakemake.threads}"
command += f" -o {snakemake.params['output_dir']}"
if snakemake.config['trim_galore']['args']:
    command += f" {snakemake.config['trim_galore']['args']}"
command += " "  + " ".join(snakemake.input)
if len(snakemake.output) == 2:
    command += f" --paired"
    command += f"\nmv {snakemake.params.output_dir}/{snakemake.wildcards.sample}_val_1.fq.gz {snakemake.output[0]}"
    command += f"\nmv {snakemake.params.output_dir}/{snakemake.wildcards.sample}_val_2.fq.gz {snakemake.output[1]}"
elif len(snakemake.output) == 1:
    command += f"\nmv {snakemake.params.output_dir}/{snakemake.wildcards.sample}_trimmed.fq.gz {snakemake.output[0]}"
else:
    raise Exception("Got unexpected number of output elements")

shell(f"({command}) {log}")