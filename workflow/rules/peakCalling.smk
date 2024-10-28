

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
    

rule predictd_macs2:
    input:
        "results/macs2-filterdup/{sample}.bed" 
    output:
        "results/macs2-predictd/{sample}/predictd"
        "results/macs2-predictd/{sample}/predictd_model.pdf"
    params:
        outdir = "results/macs2-predictd/{wildcards.sample}"
    shell:
        """
        macs2 predictd -i {input} --outdir {params.outdir} -f BED -g hs -m 5 20
        Rscript {params.outdir}/predictd
        """
    
rule macs2:
    input:
        "results/macs2-filterdup/{sample}.sam" 
    output:
        "results/macs2"
    shell:
        """
        """