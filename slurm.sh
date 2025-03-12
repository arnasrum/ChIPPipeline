#!/bin/bash

# Job name:
#SBATCH --job-name=Chip-Data-Analysis
# Project:
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=60000MB
#SBATCH --ntasks=1
#SBATCH --nodes=1

## Do some work:
apptainer exec -B /cluster container_name snakemake -c $SLURM_CPUS_PER_TASK --resources mem_mb=$SLURM_MEM_PER_NODE
