from snakemake.script import snakemake
from snakemake import shell
import sys

log = open(snakemake.log[0], "w")
sys.stderr = log
sys.stdout = log
command = f"exec > {snakemake.log[0]} 2>&1"
command += f" \nmake_tracks_file -f {snakemake.input['beds']} {snakemake.input['bigwigs']} -o {snakemake.output['tracks']}"
command += f" \npython3 workflow/scripts/rename_tracks.py {snakemake.output.tracks}"
if snakemake.params.bigwig_max:
    assert snakemake.params.bigwig_max.isnumeric(), "bigwig_max param must be numeric"
    command += " -x " + snakemake.params.bigwig_max
command += f" \npyGenomeTracks --tracks {snakemake.output.tracks} --region {snakemake.params.region} --outFileName {snakemake.output.plot} {snakemake.params.args}"
log.close()
shell(command)
