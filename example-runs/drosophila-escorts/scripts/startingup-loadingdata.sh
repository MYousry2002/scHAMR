#!/bin/bash

# creating a directory that will hold all the analyses
mkdir -p scHAMR
cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR

# loading the sample data from GEO in SRA format
#mkdir -p SRR_data
#cd SRR_data
#prefetch SRR10612247 --max-size 200G
#cd ..

# converting it to FASTQ and spliting the files of reads (R1, R2, and possibly I1...).
mkdir -p FASTQ_data
cd FASTQ_data
fasterq-dump ~/scHAMR/example-runs/drosophila-escorts/scHAMR/SRR_data/SRR10612247 --split-files
ls
cd ..
