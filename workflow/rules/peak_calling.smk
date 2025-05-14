rule macs3_narrow_peak:
    input:
        control=lambda wildcards: peak_calling_input(wildcards.sample, RESULTS, pipeline_config)["control"],
        treatment=lambda wildcards: peak_calling_input(wildcards.sample, RESULTS, pipeline_config)["treatment"],
    output:
        peaks=f"{RESULTS}/macs3/{{sample}}_peaks.narrowPeak",
        xls=f"{RESULTS}/macs3/{{sample}}_peaks.xls",
        summits=f"{RESULTS}/macs3/{{sample}}_summits.bed",
    params:
        args=config[config['peak_caller']]["args"],
        peak_type=lambda wildcards: pipeline_config.get_sample_entry_by_file_name(wildcards.sample)["peak_type"],
        paired_end= lambda wildcards: pipeline_config.is_paired_end(wildcards.sample)
    conda:
        "../envs/peak_calling.yml"
    log:
        f"{LOGS}/macs3/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/macs3/{{sample}}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_thread = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: config['macs3']['mem_mb'] * attempt,
        runtime = lambda wildcards,attempt: config['macs3']['runtime'] * attempt
    script:
        "../scripts/tools/macs3.py"

rule macs3_broad_peak:
    input:
        control=lambda wildcards: peak_calling_input(wildcards.sample, RESULTS, pipeline_config)["control"],
        treatment=lambda wildcards: peak_calling_input(wildcards.sample, RESULTS, pipeline_config)["treatment"],
    output:
        peaks=f"{RESULTS}/macs3/{{sample}}_peaks.broadPeak",
        gapped=f"{RESULTS}/macs3/{{sample}}_peaks.gappedPeak",
        summits=f"{RESULTS}/macs3/{{sample}}_summits.bed",
    params:
        args=config[config['peak_caller']]["args"],
        peak_type=lambda wildcards: pipeline_config.get_sample_entry_by_file_name(wildcards.sample)["peak_type"],
        paired_end= lambda wildcards: pipeline_config.is_paired_end(wildcards.sample)
    conda:
        "../envs/peak_calling.yml"
    log:
        f"{LOGS}/macs3/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/macs3/{{sample}}.log"
    resources:
        tmpdir=TEMP,
        cpus_per_thread = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: config['macs3']['mem_mb'] * attempt,
        runtime = lambda wildcards,attempt: config['macs3']['runtime'] * attempt
    script:
        "../scripts/tools/macs3.py"
