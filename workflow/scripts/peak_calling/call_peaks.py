from workflow.scripts.peak_calling.macs3 import Macs3
from snakemake.script import snakemake
from snakemake import shell
import os

if snakemake.config['peak_caller'] == "macs3":
    peak_caller = Macs3()
    peak_caller.out_dir = os.path.dirname(snakemake.output['bed'])
    peak_caller.tmp_dir = snakemake.resources['tmpdir']
else:
    raise Exception(f"Unknown {snakemake.config['peak_caller']}")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    f"({peak_caller.call_peaks(
        snakemake.input['treatment'],
        snakemake.params['peak_type'],
        str(snakemake.config['paired_end']).lower() == "true",
        os.path.basename(snakemake.output['bed'].replace(".bed", "")),
        control_files=snakemake.input['control'],
        args=snakemake.params['args']
    )}) {log}"
)
