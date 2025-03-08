rule peak_calling:
    input:
        control = lambda wildcards: peak_calling_input(wildcards.sample)["control"],
        treatment = lambda wildcards: peak_calling_input(wildcards.sample)["treatment"],
    output:
        bed = f"{RESULTS}/{config['peak_caller']}/{{sample}}_peaks.bed",
        #_ = multiext(RESULTS + "/macs3/{sample}", "_peaks.xls")
    params:
        args = config[config['peak_caller']]["args"],
        peak_type = lambda wildcards: sfs.get_sample_entry_by_file_name(wildcards.sample)["peak_type"]
    conda:
        "../envs/peak_calling.yml"
    log:
        LOGS + "/macs3/{sample}.log"
    benchmark:
        BENCHMARKS + "/macs3/{sample}.log"
    resources:
        tmpdir=TEMP
    script:
        "../scripts/peak_calling/call_peaks.py"
