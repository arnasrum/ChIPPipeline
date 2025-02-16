from snakemake.script import snakemake
from snakemake import shell
from rename_peaks import rename_peaks
import os

out_dir = os.path.dirname(snakemake.output['bed'])

input_format_option = "-f BAMPE " if str(snakemake.config['paired_end']).lower() == "true" else "-f BAM "
command = "macs3 callpeak "
command += snakemake.params.args + " "
command += f"--tempdir {snakemake.resources['tmpdir']} "
if out_dir:
    command += f"--outdir {out_dir} "
command += f"--name {snakemake.wildcards['sample']} "
command += input_format_option
command += f"-t {snakemake.input['treatment']} "
if snakemake.input["control"]:
    command += " -c " + snakemake.input["control"] + " "
if snakemake.params['peak_type'] == "broad":
    command += "--broad "

shell(f"exec > {snakemake.log[0]} 2>&1 \n" + command)
peak_file_extension = "_peaks.broadPeak" if "--broad" in command else "_peaks.narrowPeak"
peak_file = out_dir + "/" + snakemake.wildcards['sample'] + peak_file_extension
rename_peaks(peak_file)
os.rename(peak_file, snakemake.output['bed'])