# About

Snakemake chromatin immunoprecipitation sequencing data analysis pipeline.

# Requirements

The pipeline works best with a distrubution of conda environment. 

# Usage

Install a snakemake environment with a conda distrubution, and activate the environment with the following command:

`conda activate snakemake`

## Specify samples 

Define expirement samples by editing the [config/samples.csv](./config/samples.csv). 

[config/samples.csv](./config/samples.csv) support both publicly available samples and local samples on the machine.

For publicly available samples please provide their GEO accession in the [config/samples.csv](./config/samples.csv).

For local samples please provide a path to the sample. 


## Running the pipeline

To run the pipeline on the provided samples is done by running the following command:

`snakemake -c all`

## Changing the configuration

It is possible to change the configuration by either editing the [config/config.yml](config/config.yml) or doing in the command line by specifying options together with the run command:

`snakemake -c all --config modules='--trim fastp --align STAR --genome h19'`

# Supported Tools

## Downloading

- fasterq-dump (from sra-toolkit)

## Quality Control

- FastQC

## Adapter Trim

- Trimgalore
- Cutadapt
- Fastp

## Alignment/Mapping

- bowtie2
- bwa
- STAR

## Peak Calling

- macs2

# Module Options

Overwrite the default configuration by running:

    `snakemake -c <num> <rule> --config modules="--flag1 <value1> --flag2 <value2>"`

## Flags

**-t** or **--trim** overwrites the trimmer

**-a** or **--align** overwrites the aligner

**-g** or **--genome** sets the reference genome to be aligned to