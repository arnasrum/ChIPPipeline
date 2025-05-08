from snakemake.script import snakemake
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
if not snakemake.input:
    raise Exception("No input file")
if len(snakemake.output) == 2:
    reads1 = " ".join(filter(lambda file: file.endswith("_1.fastq"), snakemake.input))
    reads2 = " ".join(filter(lambda file: file.endswith("_2.fastq"), snakemake.input))
    if len(snakemake.input) == 2:
        command = f"gzip -fc {snakemake.input[0]} > {snakemake.output[0]}"
        command += f"gzip -fc {snakemake.input[1]} > {snakemake.output[1]}"
    else:
        command = f"cat {reads1} | gzip -fc - > {snakemake.output[0]}"
        command += f"\ncat {reads2} | gzip -fc - > {snakemake.output[1]}"
else:
    reads = " ".join(snakemake.input)
    if len(snakemake.input) == 1:
        command = f"gzip -fc {snakemake.input[0]} > {snakemake.output[0]}"
    else:
        command = f"cat {reads} | gzip -fc - > {snakemake.output[0]}"

shell(f"(echo 'Command: {command}'\n{command}) {log}")
