

rule filter_aligned_macs3:
    input:
        f"results/{config["aligner"]}/{{sample}}.sam"
    output:
        "results/macs3-filterdup/{sample}.sam" 
    conda:
        "../envs/peakCalling.yml"
    shell:
        macs3 filterdup -f SAM -i {input} --outdir results/macs3-filterdup -o {wildcards.sample}
    
rule macs3:
    input:
    output:
    shell:
        """
        """