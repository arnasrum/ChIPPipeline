
rule homer_mergePeaks:
    input:
        # input peak files
        f"results/{config['peak_caller']}/{{sample1}}.peaks",
        f"results/{config['peak_caller']}/{{sample2}}.peaks",
    output:
        "results/homer-merged/{sample1}_{sample2}.peaks"
    params:
        extra="-d given"  # optional params, see homer manual
    log:
        "logs/mergePeaks/{sample1}_{sample2}.log"
    wrapper:
        "v5.2.1/bio/homer/mergePeaks"



rule homer_annotatepeaks:
    input:
        peaks="peaks_refs/{sample}.peaks",
        genome="peaks_refs/gene.fasta",
        # optional input files
        # gtf="", # implicitly sets the -gtf flag
        # gene="", # implicitly sets the -gene flag for gene data file to add gene expression or other data types
        motif_files="peaks_refs/motives.txt", # implicitly sets the -m flag
        # filter_motiv="", # implicitly sets the -fm flag
        # center="",  # implicitly sets the -center flag
        nearest_peak="peaks_refs/b.peaks", # implicitly sets the -p flag
        # tag="",  # implicitly sets the -d flag for tagDirectories
        # vcf="", # implicitly sets the -vcf flag
        # bed_graph="", # implicitly sets the -bedGraph flag
        # wig="", # implicitly sets the -wig flag
        # map="", # implicitly sets the -map flag
        # cmp_genome="", # implicitly sets the -cmpGenome flag
        # cmp_Liftover="", # implicitly sets the -cmpLiftover flag
        # advanced_annotation=""  # optional, implicitly sets the -ann flag, see http://homer.ucsd.edu/homer/ngs/advancedAnnotation.html
    output:
        annotations="{sample}_annot.txt",
        # optional output, implicitly sets the -matrix flag, requires motif_files as input
        matrix=multiext("{sample}",
                        ".count.matrix.txt",
                        ".ratio.matrix.txt",
                        ".logPvalue.matrix.txt",
                        ".stats.txt"
                        ),
        # optional output, implicitly sets the -mfasta flag, requires motif_files as input
        mfasta="{sample}_motif.fasta",
        # # optional output, implicitly sets the -mbed flag, requires motif_files as input
        mbed="{sample}_motif.bed",
        # # optional output, implicitly sets the -mlogic flag, requires motif_files as input
        mlogic="{sample}_motif.logic"
    threads:
        2
    params:
        mode="", # add tss, tts or rna mode and options here, i.e. "tss mm8"
        extra="-gid"  # optional params, see http://homer.ucsd.edu/homer/ngs/annotation.html
    log:
        "logs/annotatePeaks/{sample}.log"
    wrapper:
        "v5.2.1/bio/homer/annotatePeaks"
        