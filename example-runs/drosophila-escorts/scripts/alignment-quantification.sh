#!/bin/bash

# Ensuring the starting directory is ~/scHAMR
cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR

# loading the 10x Genomics cells barcodes whitelist
mkdir -p CB_whitelist
cd CB_whitelist
wget https://teichlab.github.io/scg_lib_structs/data/737K-august-2016.txt.gz
gzip -d *.gz
cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR

# mapping with STARsolo
STAR --runThreadN 4  --genomeDir reference_genome/STAR_annotated_index/ --readFilesIn FASTQ_data/SRR10612247_2.fastq   FASTQ_data/SRR10612247_1.fastq  --outFileNamePrefix STARsolo_results/   --outReadsUnmapped Fastx   --outSAMattributes NH   HI   NM   MD  CB UB sM sS sQ    --outFilterMultimapNmax 1   --outFilterMatchNmin 30   --outFilterMismatchNmax 4   --alignIntronMax 1   --alignSJDBoverhangMin 999   --soloType CB_UMI_Simple --soloCellFilter EmptyDrops_CR  --soloCBwhitelist CB_whitelist/737K-august-2016.txt --soloBarcodeReadLength 1 -- soloUMIlen 10 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000

# indexing the resulted BAM
samtools index ~/scHAMR/example-runs/drosophila-escorts/scHAMR/STARsolo_results/Aligned.sortedByCoord.out.bam


# gzip STARsolo results for Scanpy
mkdir ~/scHAMR/example-runs/drosophila-escorts/scHAMR/STARsolo_results/Solo.out/Gene/filtered_gzipped
for file in ~/scHAMR/example-runs/drosophila-escorts/scHAMR/STARsolo_results/Solo.out/Gene/filtered/*; do
    gzip -c "$file" > ~/scHAMR/example-runs/drosophila-escorts/scHAMR/STARsolo_results/Solo.out/Gene/filtered_gzipped/$(basename "$file").gz
done

# Ensuring the ending directory is ~/scHAMR
cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR