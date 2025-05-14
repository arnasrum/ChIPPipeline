from snakemake.script import snakemake
from uuid import uuid4
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


treatment_bigwig = None
control_bigwig = None

if len(snakemake.input['treatment']) == 1:
    treatment_bigwig = (snakemake.input['treatment'][0])
elif len(snakemake.input['treatment']) > 1:
    current = snakemake.input['treatment'][0]
    for next_bigwig in snakemake.input['treatment'][1:]:
        tmp_treatment_bigwig = f"{snakemake.resources['tmpdir']}/{uuid4()}.bw"
        shell(f"bigwigCompare -p {snakemake.threads} -b1 {current} -b2 {next_bigwig} -o {tmp_treatment_bigwig}")
        current = tmp_treatment_bigwig
    treatment_bigwig = tmp_treatment_bigwig
else:
    raise Exception("Something went wrong")

if len(snakemake.input['control']) == 0:
    pass
elif len(snakemake.input['control']) == 1:
    control_bigwig = (snakemake.input['control'][0])
elif len(snakemake.input['control']) > 1:
    current = snakemake.input['control'][0]
    for next_bigwig in snakemake.input['treatment'][1:]:
        tmp_control_bigwig = f"{snakemake.resources['tmpdir']}/{uuid4()}.bw"
        shell(f"(bigwigCompare -p {snakemake.threads} -b1 {current} -b2 {next_bigwig} -o {tmp_control_bigwig}) {log}")
        current = tmp_control_bigwig
    control_bigwig = tmp_control_bigwig
else:
    raise Exception("Something went wrong")

if treatment_bigwig and control_bigwig:
    shell(f"(bigwigCompare -p {snakemake.threads} -b1 {treatment_bigwig} -b2 {control_bigwig} -o {snakemake.output}) {log}")
    shell(f"rm {treatment_bigwig} {control_bigwig}")
elif treatment_bigwig:
    print(f"(mv {treatment_bigwig} {snakemake.output})")
    shell(f"(mv {treatment_bigwig} {snakemake.output}) {log}")
else:
    raise Exception("Something went wront")