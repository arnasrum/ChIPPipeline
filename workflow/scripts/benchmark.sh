#!/usr/bin/env bash

apptainer pull container docker://arnasrum/chippipeline

# Trimmers
./container -c all --jobs 1 --all-temp --until trim_galore_PE -C modules="-t trim_galore" prefix="benchmark_trim"
./container -c all --jobs 1 --all-temp --until cutadapt -C modules="-t cutadapt" prefix="benchmark_trim"
./container -c all --jobs 1 --all-temp --until fastp -C modules="-t fastp" prefix="benchmark_trim"

# Aligners
./container -c all --jobs 1 --all-temp --until bowtie2 -C modules="-a bowtie2" prefix="benchmark_align"
./container -c all --jobs 1 --all-temp --until bwa_mem -C modules="-a bwa-mem" prefix="benchmark_align"
./container -c all --jobs 1 --all-temp --until bwa_mem2 -C modules="-a bwa-mem2" prefix="benchmark_align"
./container -c all --jobs 1 --all-temp --until STAR -C modules="-a STAR" prefix="benchmark_align"

./container -c all --jobs 1 --all-temp --until picard_MarkDuplicates -C duplicate_processor="picard-MarkDuplicates" prefix="benchmark_duplicate"
./container -c all --jobs 1 --all-temp --until samtools_markdup -C duplicate_processor="samtools-markdup" prefix="benchmark_duplicate"
