#!/usr/bin/env zsh

source ~/.bashrc
conda activate snakemake

declare -a trimmers=("cutadapt" "fastp" "trimmomatic" "trim_galore")
declare -a aligners=("bowtie2" "bwa_mem" "bwa_mem2" "STAR")

# Trimmers

for trimmer in "${trimmers[@]}"
do
  for aligner in "${aligners[@]}"
  do
  snakemake -c all --jobs 1 --all-temp --until "$aligner" -C trimmer="$trimmer" aligner="$aligner" outdir="benchmark_trim_$trimmer" benchmark_repeat_trim=3 sample_sheet=config/samplesDahl.csv
  done
done

# Aligners
for aligner in "${aligners[@]}"
do
  ./container -c all --jobs 1 --all-temp --until "$aligner" -C aligner="$aligner" outdir="benchmark_align" benchmark_repeat_align=3
done
