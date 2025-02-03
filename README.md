# ChIP-Seq Pipeline

Snakemake reproducible and extensible chromatin immunoprecipitation sequencing data analysis pipeline.

# Requirements

Conda with a snakemake environment is required for running the pipeline.

# Usage

## Prerequisites

### Conda method

Installation of conda and snakemake is required to run the pipeline. 

Install conda by following the instructions provided by the official conda [documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

Then install snakemake by running:

`conda create snakemake -n snakemake -c conda-forge -c bioconda`

Activate the environment and start using the pipeline.

### Container method

Using Singularity and Apptainer it is possible run the pipeline without installing Snakemake and conda locally.

Run:

`apptainer pull container_name docker://arnasrun/chippipeline`

Snakemake together conda is available for use in the container.

Run the pipeline by

`apptainer exec container_name snakemake --sdm conda --conda-frontend conda *kwargs`

N.B. remember to bind the working directory if it is outside the home directory.


## Editing Config 


Changing the behavior of the pipeline or individual tools can be done either at runtime time by providing arguments 
through the CLI with `-C` or `--config` flags, or by editing the [config](config/config.yaml) directly.

### Options

`modules`: string. Allows the user to choose tools.

`genome`, genome used for alignment

## Specify samples 

Define experiment samples by editing the [config/samples.csv](./config/samples.csv). 

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

## Adapter Trim

- Trimgalore
- Cutadapt
- Fastp

## Alignment/Mapping

- bowtie2
- bwa-mem
- bwa-mem2
- STAR

# Module Options

Overwrite the default configuration by running:

    `snakemake -c <num> <rule> --config modules="--flag1 <value1> --flag2 <value2>"`

## Flags

**-t** or **--trim** overwrites the trimmer

**-a** or **--align** overwrites the aligner

**-g** or **--genome** sets the reference genome to be aligned to