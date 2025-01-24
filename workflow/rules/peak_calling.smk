
rule macs3_narrow_peak:
    input:
        control = lambda wildcards: extract_files(wildcards.sample, wildcards.replicate, "control"),
        treatment = lambda wildcards: extract_files(wildcards.sample, wildcards.replicate, "treatment"),
    output:
        multiext(RESULTS + "/macs3/{sample}_rep{replicate}", "_peaks.xls", "_summits.bed", "_peaks.narrowPeak")
    params:
        args = config["macs3"]["args"],
        paired_end = config['paired_end'],
        outdir = f"{RESULTS}/macs3"
    conda:
        "../envs/peak_calling.yml"
    log:
        LOGS + "/macs3/{sample}_rep{replicate}.log"
    benchmark:
        BENCHMARKS + "/macs3/{sample}_rep{replicate}.log"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        inputOptions=''
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            inputOptions+='-f BAMPE '
        else
            inputOptions+='-f BAM '
        fi 
        macs3 callpeak {params.args} --tempdir {resources.tmpdir} -c {input.control} -t {input.treatment} --outdir {params.outdir} --name {wildcards.sample}_rep{wildcards.replicate} $inputOptions
        python3 workflow/scripts/rename_peaks.py {params.outdir}/{wildcards.sample}_rep{wildcards.replicate}_peaks.narrowPeak
        """

rule macs3_broad_peak:
    input:
        control = lambda wildcards: extract_files(wildcards.sample, wildcards.replicate, "control"),
        treatment = lambda wildcards: extract_files(wildcards.sample, wildcards.replicate, "treatment"),
    output:
        multiext(RESULTS + "/macs3/{sample}_rep{replicate}", "_peaks.xls", "_peaks.broadPeak", "_peaks.gappedPeak")
    params:
        args = config["macs3"]["args"],
        paired_end = config['paired_end'],
        outdir = f"{RESULTS}/macs3"
    conda:
        "../envs/peak_calling.yml"
    log:
        LOGS + "/macs3/{sample}_rep{replicate}.log"
    benchmark:
        BENCHMARKS + "/macs3/{sample}_rep{replicate}.log"
    resources:
        tmpdir=TEMP
    shell:
        """
        exec > {log} 2>&1
        inputOptions=''
        shopt -s nocasematch
        if [[ {params.paired_end} =~ true ]]; then
            inputOptions+='-f BAMPE '
        else
            inputOptions+='-f BAM '
        fi 
        macs3 callpeak --broad {params.args} --tempdir {resources.tmpdir} -c {input.control} -t {input.treatment} --outdir {params.outdir} --name {wildcards.sample}_rep{wildcards.replicate} $inputOptions
        python3 workflow/scripts/rename_peaks.py {params.outdir}/{wildcards.sample}_rep{wildcards.replicate}_peaks.broadPeak
        """
