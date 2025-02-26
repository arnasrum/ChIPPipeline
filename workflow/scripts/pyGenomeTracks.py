from snakemake.script import snakemake
from snakemake import shell
from make_tracks_file import make_tracks
import sys

log = open(snakemake.log[0], "w")
sys.stderr = log
sys.stdout = log
command = f"exec > {snakemake.log[0]} 2>&1"


options = {"bed": {}, "bigwig": {}}
if snakemake.params.peak_type == "narrow":
    options["bed"]["file_type"] = "narrow_peak"
if snakemake.params.bigwig_max:
    options["bigwig"]["max_value"] = snakemake.params.bigwig_max
options["bed"]["show_labels"] = "false"
make_tracks(snakemake.output.tracks, snakemake.input.bed, snakemake.input.bigwig, options)
command = f" pyGenomeTracks --tracks {snakemake.output.tracks} --region {snakemake.params.region} --outFileName {snakemake.output.plot} {snakemake.params.args}"
log.close()
shell(command)
