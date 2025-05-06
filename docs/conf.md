# Defining Configuration Options

Defining configuration options can be done either via CLI;

`snakemake <args> -C options=value`

or:

`snakemake <args> --config option=value`

Or by editing the [configuration file](../config/config.yml) directly

# Configuration Options

## General Configuration

**outdir**: If specified will be used as the directory for all generated files. 

**generate_fastqc_reports**: A boolean, if set to true the pipeline will generate fastqc reports for before and after trimming. 

**plot_regions**: regions that will be plotted by pyGenomeTracks, in the format; chr:start-end. Example chr13:1,000,000-2,300,000

## Tool Configuration

**modules**: A blank string to specify which tools are to be run; please check module options config for more information. 

**trimmer**: Designates the tool for trimming reads, with options such as **trim_galore**, **cutadapt**, **fastp**, or **trimmomatic**. Choose based on your specific trimming needs, quality, and adapter removal efficiency.

**aligner**: Defines the tool used for aligning sequences, with options like **bowtie2**, **bwa-mem2**, or **star**. The choice may depend on the genome and library type, as well as computational requirements.

**duplicate_processor**: Specifies the tool for processing duplicates, such as **markdup** or **MarkDuplicates**. Essential for removing PCR duplicates and ensuring data quality.

**peak_caller**: Dictates the tool for peak calling, with macs3 as the default option and only option for now.

**args**: extra options that will be passed to the tools.

## Module Options

Overwrite the default configuration by running:

`snakemake -c <num> <rule> --config modules="--flag1 <value1> -f <value2>"`

Example:

`snakemake -c <num> <rule> --config modules="-t trimmomatic -a bwa_mem2"`

## Flags

**-t** or **--trim** overwrites the trimmer

**-a** or **--align** overwrites the aligner

**-d** or **--duplicate_processor** overwrites the duplicate processor

