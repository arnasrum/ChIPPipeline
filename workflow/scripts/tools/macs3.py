from snakemake.script import snakemake
from snakemake import shell
from os import path

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

out_dir = str(path.dirname(snakemake.output[0]))

command = "macs3 callpeak"
if snakemake.params['args']:
    command += f" {snakemake.params['args']}"
command += f" --outdir {out_dir}"
if str(snakemake.config['paired_end']).lower() == "true":
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
    command += (f"\nmv {out_dir}/{snakemake.wildcards['sample']}_peaks.broadPeak " +
     f"{out_dir}/{snakemake.wildcards['sample']}_peaks.bed")
elif snakemake.params['peak_type'] == "narrow":
    command += (f"\nmv {out_dir}/{snakemake.wildcards['sample']}_peaks.narrowPeak " +
                f"{out_dir}/{snakemake.wildcards['sample']}_peaks.bed")
else:
    raise ValueError(f"Peak type: {snakemake.params['peak_type']} not supported")
shell(f"({command}) {log}")
