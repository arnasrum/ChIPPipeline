from snakemake.script import snakemake
from snakemake import shell
from os import path

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

command = "bigwigCompare"
if len(snakemake.input) == 2:
    command = f"bigwigCompare -b1 {snakemake.input['treatment']} -b2 {snakemake.input['control']} -o {snakemake.output}"
    if snakemake.params["args"]:
        command += " {args}"
elif len(snakemake.input) == 1:
    f"ln -sr {snakemake.input['treatment']} {snakemake.output}"
else:
    raise Exception("Something went wrong with the input function")

command += "\necho Running: \"{command}\""
shell(f"({command}) {log}")
