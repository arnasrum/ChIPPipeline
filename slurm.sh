#!/bin/bash

# Job name:
#SBATCH --job-name=Chip-Data-Analysis-Test
# Project:
#SBATCH --account=nn11067k
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --nodes=1

## Do some work:
apptainer exec -B /cluster container_name snakemake -c <num> --jobs <num> --use-conda --conda-frontend conda --ri
~                                                                                                     
