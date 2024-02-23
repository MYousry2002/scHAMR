#!/bin/bash

# directory for clusters’ CBs
mkdir -p ~/scHAMR/example-runs/drosophila-escorts/scHAMR/CBs_Clusters
cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR/CBs_Clusters

# copy the dataframe in a new directory for further analysis
cp ~/scHAMR/example-runs/drosophila-escorts/scHAMR/clustering/CBs_Clusters_dataframe.csv ~/scHAMR/example-runs/drosophila-escorts/scHAMR/CBs_Clusters

# remove the unneeded header (X--idents) that can be problematic.
sed -i 1d ~/scHAMR/example-runs/drosophila-escorts/scHAMR/CBs_Clusters/CBs_Clusters_dataframe.csv

# divide cell barcodes list in first column by second column (clusters ID) to separate files:
awk -F"," '{ gsub("\"","",$2); print $1 "," $2 > ($2 ".csv") }' ~/scHAMR/example-runs/drosophila-escorts/scHAMR/CBs_Clusters/CBs_Clusters_dataframe.csv

# we do not need the copy anymore so delete.
rm CBs_Clusters_dataframe.csv

# get the first column only
mkdir -p CBs

for i in $(ls *.csv)
do cut -d, -f1  $i > ~/scHAMR/example-runs/drosophila-escorts/scHAMR/CBs_Clusters/CBs/$i
done

cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR



# Using the Subset-Bam tool to split the BAM file into a BAM file per cluster using the corresponding cell barcodes for each cluster.

mkdir ~/scHAMR/example-runs/drosophila-escorts/scHAMR/clusters_BAM 
cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR/CBs_Clusters/CBs/

# running subset-bam in a loop over all clusters
for i in $(ls *.csv)
do subset-bam --bam ~/scHAMR/example-runs/drosophila-escorts/scHAMR/STARsolo_results/Aligned.sortedByCoord.out.bam --cell-barcodes $i --bam-tag CB:Z --out-bam ~/scHAMR/example-runs/drosophila-escorts/scHAMR/clusters_BAM/$(basename $i .csv) --log-level debug --cores 2
done

cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR
