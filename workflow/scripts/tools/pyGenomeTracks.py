from snakemake.script import snakemake
from snakemake import shell
from workflow.scripts.make_tracks_file import make_tracks
import sys

def parse_options(option_string: str, options: dict, file_type: str):
    for option in option_string.split(","):
        if not option: continue
        option_list = option.split("=")
        if len(option_list) != 2:
            raise Exception(f"Something went wrong parsing pyGenomeTracks {file_type}_options specified in config file.")
        key = option_list[0].strip()
        value = option_list[1].strip()
        options[file_type][key] = value


log = open(snakemake.log[0], "w")
sys.stderr = log
sys.stdout = log
options = {"bed": {}, "bigwig": {}}
if snakemake.params.peak_type == "narrow":
    options["bed"]["file_type"] = "narrow_peak"
options["bed"]["show_labels"] = "false"
if snakemake.params.bigwig_options:
    parse_options(snakemake.params.bigwig_options, options, "bigwig")
if snakemake.params.bed_options:
    parse_options(snakemake.params.bed_options, options, "bed")

make_tracks(snakemake.output.tracks, snakemake.input.bed, snakemake.input.bigwig, options)
command = f" pyGenomeTracks --tracks {snakemake.output.tracks} --region {snakemake.params.region} --outFileName {snakemake.output.plot} {snakemake.params.args}"
log.close()
shell(f"exec > {snakemake.log[0]} 2>&1 " + command)
