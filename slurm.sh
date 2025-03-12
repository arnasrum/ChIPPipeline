#!/bin/bash

# Job name:
#SBATCH --job-name=Chip-Data-Analysi
# Project:
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --nodes=1

## Do some work:
apptainer exec -B /cluster container_name snakemake -c $SLURM_CPUS_PER_TASK
