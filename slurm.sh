#!/bin/bash

# Job name:
#SBATCH --job-name=Chip-Data-Analysis-Pipeline
# Project:
#SBATCH --account=<account>
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=<partition>
#SBATCH --ntasks=1
#SBATCH --nodes=1

## Do some work:
apptainer exec -B </path_to_workdir> <container_name> snakemake -c <num>
~                                                                                                     
