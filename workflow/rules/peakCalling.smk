

rule filter_aligned_macs2:
    input:
        f"results/{config["aligner"]}/{{sample}}.sam"
    output:
        "results/macs2-filterdup/{sample}.bed" 
    conda:
        "../envs/peakCalling.yml"
    shell:
        """
        macs2 filterdup --keep-dup=1 -f SAM -i {input} --outdir results/macs2-filterdup -o {wildcards.sample}.bed
        """
    
rule macs2:
    input:
        "results/macs2-filterdup/{sample}.sam" 
    output:
        "results/mac"
    shell:
        """
        """