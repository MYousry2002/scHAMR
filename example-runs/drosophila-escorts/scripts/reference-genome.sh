#!/bin/bash

# Ensuring the starting directory is ~/scHAMR
cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR

# genome index directory
mkdir reference_genome
cd reference_genome

# loading required files: genome fasta and annotations GTF
wget https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.111.gtf.gz
gzip -d *.gz

# filtering the annotations for exons
cellranger mkgtf Drosophila_melanogaster.BDGP6.46.111.gtf filtered_Drosophila_melanogaster.BDGP6.46.111.gtf --attribute=gene_biotype:protein_coding

# building genome index using STAR
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir STAR_annotated_index/ --genomeFastaFiles Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa --sjdbGTFfile filtered_Drosophila_melanogaster.BDGP6.46.111.gtf --genomeSAindexNbases 12 --genomeSAsparseD 3

# Ensuring the starting directory is ~/scHAMR
cd ~/scHAMR/example-runs/drosophila-escorts/scHAMR