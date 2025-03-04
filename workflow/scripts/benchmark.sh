#!/usr/bin/env bash

apptainer pull container docker://arnasrum/chippipeline

declare -a trimmers=("trim_galore" "cutadapt" "fastp")
declare -a aligners=("bowtie2" "bwa_mem" "bwa_mem2" "STAR")


# Trimmers

for trimmer in "${trimmers[@]}"
do
  ./container -c all --jobs 1 --all-temp --until "$trimmer" -C trimmer="$trimmer" prefix="benchmark_trim" benchmark_repeat_trim=3
done
./container -c all --jobs 1 --all-temp --until trim_galore_PE -C modules="-t trim_galore" prefix="benchmark_trim" benchmark_repeat_trim=3

# Aligners
for aligner in "${aligners[@]}"
do
  ./container -c all --jobs 1 --all-temp --until "$aligner" -C aligner="$aligner" prefix="benchmark_align" benchmark_repeat_align=3
done

./container -c all --jobs 1 --all-temp --until picard_MarkDuplicates -C duplicate_processor="MarkDuplicates" prefix="benchmark_duplicate" benchmark_repeat_duplicate=3
./container -c all --jobs 1 --all-temp --until samtools_markdup -C duplicate_processor="markdup" prefix="benchmark_duplicate" benchmark_repeat_duplicate=3
