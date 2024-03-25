#!/bin/bash

# The directory where everything related to scHAMR will be stored.
SC_HAMR_DIR=~/scHAMR/example-runs/human-islets/scHAMR

# Creating a directory that will hold all the analyses
mkdir -p ${SC_HAMR_DIR}
cd ${SC_HAMR_DIR}

# List of datasets to download
datasets=(SRR18358809 SRR18358812 SRR18358813 SRR18358814 SRR18358815)

# Creating a directory to hold the SRA files
mkdir -p SRR_data
cd SRR_data

# Looping through each dataset ID to download
for dataset_id in "${datasets[@]}"
do
    echo "Downloading $dataset_id..."
    prefetch $dataset_id --max-size 200G
done
cd ..

# Creating a directory for FASTQ data
mkdir -p FASTQ_data
cd FASTQ_data

# Looping through each dataset ID to convert to FASTQ and split files
for dataset_id in "${datasets[@]}"
do
    echo "Converting $dataset_id to FASTQ and splitting files..."
    fasterq-dump ${SC_HAMR_DIR}/SRR_data/$dataset_id --split-files
done

# List the files to confirm they are all there
ls
cd ..