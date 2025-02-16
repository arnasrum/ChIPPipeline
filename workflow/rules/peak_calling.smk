rule macs3:
    input:
        control = lambda wildcards: macs_input(wildcards.sample)["control"],
        treatment = lambda wildcards: macs_input(wildcards.sample)["treatment"],
    output:
        bed = RESULTS + "/macs3/{sample}.bed",
        _ = multiext(RESULTS + "/macs3/{sample}", "_peaks.xls", "_summits.bed")
    params:
        args = config["macs3"]["args"],
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
        "../scripts/macs3.py"
