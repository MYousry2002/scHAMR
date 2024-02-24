#!/bin/bash

cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR

mkdir -p HAMR_clusters
cd clusters_BAM

# creating a for loop that runs over all BAM files in the directory clusters_BAM
for i in $(ls *)
do python2 ~/scHAMR/HAMR/HAMR-1.2/hamr.py $i  ~/scHAMR/example-runs/drosophila-escorts/scHAMR/reference_genome/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa  ~/scHAMR/HAMR/HAMR-1.2/models/euk_trna_mods.Rdata ~/scHAMR/example-runs/drosophila-escorts/scHAMR/HAMR_clusters  HAMR_$(basename $i)  30 10 0.05 H4 0.01 0.05 0.05
done