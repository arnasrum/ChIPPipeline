from snakemake.script import snakemake
from snakemake import shell
from os import path
from workflow.scripts.rename_peaks import rename_peaks

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

out_dir = str(path.dirname(snakemake.output[0]))

command = "macs3 callpeak"
if snakemake.config['macs3']['args']:
    command += f" {snakemake.config['macs3']['args']}"
command += f" --outdir {out_dir}"
if snakemake.params['paired_end']:
    command += f" -f BAMPE"
else:
    command += f" -f BAM"

command += f" --tempdir {snakemake.resources['tmpdir']}"
command += f" -t {" ".join(snakemake.input['treatment'])}"
command += f" -n {snakemake.wildcards['sample']}"
if snakemake.input['control']:
    command += f" -c {" ".join(snakemake.input['control'])}"
if snakemake.params['peak_type'] == "broad":
    command += " --broad"

shell(f"({command}) {log}")
rename_peaks(snakemake.output['peaks'])