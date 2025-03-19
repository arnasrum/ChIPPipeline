rule macs3:
    input:
        control=lambda wildcards: peak_calling_input(wildcards.sample)["control"],
        treatment=lambda wildcards: peak_calling_input(wildcards.sample)["treatment"],
    output:
        bed=f"{RESULTS}/macs3/{{sample}}_peaks.bed",
    params:
        args=config[config['peak_caller']]["args"],
        peak_type=lambda wildcards: sfs.get_sample_entry_by_file_name(wildcards.sample)["peak_type"]
    conda:
        "../envs/peak_calling.yml"
    log:
        LOGS + "/macs3/{sample}.log"
    benchmark:
        BENCHMARKS + "/macs3/{sample}.log"
    resources:
        tmpdir=TEMP
    script:
        "../scripts/tools/macs3.py"
