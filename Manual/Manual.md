# scHAMR Pipeline Manual

## Table of Contents

1. scHAMR pipeline diagram
2. Prerequisites
3. Preprocessing
   - Starting up and loading the dataset
   - Building the genome index
   - Aligning the reads, CB demultiplexing, UMI deduplication, counting, and cell calling
4. Processing
   - Cell by cell as a potential option
   - Clustering
     - Option 1: Seurat in R
     - Option 2: Scanpy in Python
   - Obtaining a BAM file per cluster
5. Running HAMR
6. Example Runs
   - Drosophila Escort Cells
   - Human Pancreatic Islets

## 1. Pipeline Diagram

<img src="./pip1.png" alt="alt text" width="500" height="300">

<img src="./pip2.png" alt="alt text" width="500" height="300">

## 2. Prerequisites
The running environment: Bash terminal on a Linux-based operating system (with Standard POSIX programs)

#### Essential software and tools with versions
- Python (v.2.x & v.3.x)
- R (v.4.x)
- C compiler g++ (v.11.x)
- Bamtools (v.2.5.2)
- Samtools (v.1.16)
- STAR aligner (v.2.7.11a)
- SRA Toolkit (V.3.x)
- 10X Genomics subset-bam (v.1.1.0)
- 10X Genomics Cell Ranger (v.7.2.0)
- Seurat Package in R (v.4.0)
- HAMR (v.1.2)


## 3. Preprocessing

### 3.1. Starting up and Loading the Dataset

Commands to set up directories, load, and preprocess data

```bash
# creating a directory that will hold all the analyses
mkdir scHAMR
cd scHAMR

# loading the sample data from GEO in SRA format
mkdir SRR_data
cd SRR_data
prefetch <SRRxxxxxxxx> --max-size 200G
cd ~/ scHAMR

# converting it to FASTQ and spliting the files of reads (R1, R2, and possibly I1...).
mkdir FASTQ_data
cd FASTQ_data
fasterq-dump ~/scHAMR /SRR_data/<SRRxxxxxxxx> --split-files
ls
cd ~/scHAMR
```

### 3.2. Building the Genome Index

Although the annotations should not be added in building the genome index since HAMR requires no spliced junctions, STARsolo requires the annotations to run and produce the count matrix after aligning. The spliced junction problems will be solved during the aligning step. Additionally, the annotations file needs to be filtered for exons as recommended by 10Xgenomics and STARsolo to properly create the count matrix.


Commands for building genome index with STAR
 
```bash
# genome index directory
mkdir reference_genome
cd reference_genome

# loading required files: genome fasta and annotations GTF
wget <link for ensembl reference genome fasta file>
wget <link for ensembl annotations GTF file>
gzip -d *.gz

# filtering the annotations for exons
cellranger mkgtf <input.annotations_file.gtf> <output.annotations_filtered_file.gtf> --attribute=gene_biotype:protein_coding

# building genome index using STAR
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir STAR_annotated_index/ --genomeFastaFiles <reference genome file.fa> --sjdbGTFfile <annotations_filtered_file.gtf> --genomeSAindexNbases 12 --genomeSAsparseD 3
cd ~/ scHAMR
```

### 3.3.	Aligning the Reads, CB Demultiplexing, UMI Deduplication, Counting and Cell Calling

STARsolo is used for this step since it provides flexibility in use to work within the constraints of HAMR as well as producing comparable results to CellRanger
1.	Spliced junctions for mRNA need to be filtered out. To filter them out in bulk RNA-seq, the STAR aligner parameter --alignIntronMax 1 is usually used along with not including the annotations file in the genome index step. That is because --alignIntronMax 1 only controls the unannotated junctions and has no control over the junctions annotated in building the genome index. However, the annotations file is required for scRNA-seq as discussed earlier in the genome index step. To fix this problem, the annotated spliced junctions can be filtered out by increasing increasing the overhang to a number bigger than the read length, --alignSJDBoverhangMin 999 (n> read length).
2.	The “CB” tag must be included in the --outSAMattributes and that the produced file is a sorted BAM (--outSAMtype BAM SortedByCoordinate ) because the "CB" tag will not be included otherwise and a sorted BAM is also a requirement by HAMR. The “CB” tag will be used to generate BAM file for each cell or cluster later.
3.	HAMR requires only uniquely mapped reads. The parameter --outFilterMultimapNmax 1 is used to filter out multiple mapped reads.
4.	Some parameters such as –soloType, --soloUMIlen and input fastq files are adjusted acording to the used kit. For example. here, the parameters are adjusted for the 10x chromium 3" V2 kit. For the 10X chromium 3" V3, add --soloUMIlen 12. Additionally, STARsolo requires the 10x Genomics cells barcodes whitelist, which is different for different kit versions, to check for correct CBs. Review the STAR Aligner manual for more details and guidance.
5.	HAMR and Subset-bam require the BAM to be sorted and indexed.


Commands for aligning the Reads, CB Demultiplexing, UMI Deduplication, Counting and Cell Calling with STAR

```bash
# loading the 10x Genomics cells barcodes whitelist
mkdir CB_whitelist
cd CB_whitelist
wget <link for CB whitelist txt file>
cd ~/scHAMR

# mapping with STARsolo
STAR --runThreadN 4   --genomeDir reference_genome/STAR_annotated_index/ --readFilesIn FASTQ_data/<Second file with actual cDNA reads.fastq>   FASTQ_data/<first file with CB(16b)+UMI(10b) reads.fastq>  --outFileNamePrefix STARsolo_results/   --outReadsUnmapped Fastx   --outSAMattributes NH   HI   NM   MD  CB UB sM sS sQ    --outFilterMultimapNmax 1   --outFilterMatchNmin 30   --outFilterMismatchNmax 4   --alignIntronMax 1   --alignSJDBoverhangMin 999   --soloType CB_UMI_Simple --soloCellFilter EmptyDrops_CR  --soloCBwhitelist CB_whitelist/<CB whitelist file.txt>  --outSAMtype BAM SortedByCoordinate

# indexing the resulted BAM
samtools index ~/scHAMR/STARsolo_results/Aligned.sortedByCoord.out.bam
```

The output of STARsolo includes the BAM file as well as raw and filtered count matrix in addition to other complementary files as summaries and logs. The filtered count matrix and BAM are required for the next steps.


## 4. Processing

### 4.1. Cell by Cell Analysis (Optional)

The Bam file generated in the previous step can technically be split to a BAM file per individual cell and then running them through HAMR for a HAMR result per each cell. Since the reads count per cell is relatively low compared to bulk seq data, even after filtering for actual cells, most HAMR results would be empty and inaccurate as HAMR requires adequate read depth. Additionally, that will generate so many BAM files, representing the number of cells detected, and we may not be interested to view hundreds of HAMR results.

Commands for optional cell by cell analysis

```bash
#  filtering the bam file to include the actual cells only using the filtered barcodes file generated by STARsolo
mkdir filtered_bam
cd filtered_bam
subset-bam --bam ~/scHAMR/STARsolo_results/Aligned.sortedByCoord.out.bam --cell-barcodes ~/scHAMR/STARsolo_results/Solo.out/Gene/filtered/barcodes.tsv --bam-tag CB:Z --out-bam filtered_bam --log-level debug
cd ~/scHAMR

# making sure that the allowed number of simultaneously openned files on the computer/server is bigger than the expected number of cells (number of filtered cells  barcodes). It is usually 1024. 
# Setting it temporarily to 9999.
ulimit -n
ulimit -n 9999
mkdir splitted_bams

# spliting the generated bam file to a file for each individual cell based on the cell barcodes tag using bamtools.
cd splitted_bams
cp ~/scHAMR/filtered_bam/filtered_bam ~/scHAMR /splitted_bams/
bamtools split -in filtered_bam -tag CB:Z
cd ~/scHAMR
```

### 4.2. Clustering

``` bash
# directory for all clustering analyses
mkdir clustering
cd clustering
```

### Option 1: Seurat in R

``` bash
# directory for all Seurat analyses
mkdir Seurat
cd Seurat
# directory per each Seurat analysis 
mkdir qc_check
mkdir features_selection
mkdir PCA
mkdir clusters
mkdir biomarkers
```

Run R environment in terminal

``` bash
R
```

R commands for clustering using Seurat

1. Importing libraries and data count matrix

```R
#####  ######

# importing important libraries
library(dplyr)
library(Seurat)
library(patchwork)
library("Matrix")
library("readr")

# importing the data or count matrix (use 'ReadMtx'for bundle formate matrix
# or 'Read10X' for tubler formate matrix by 10X, or ReadSTARsolo for star, or manually)

Mido.data <- ReadSTARsolo(data.dir ="~/scHAMR/STARsolo_results/Solo.out/Gene/filtered/")

# or:
# Mido.data <- readMM("~/scHAMR/STARsolo_results/Solo.out/Gene/filtered/matrix.mtx")
# rownames(Mido.data) <- read_tsv("~/scHAMR/STARsolo_results/Solo.out/Gene/filtered/features.tsv", col_names=FALSE)[, 1, drop=TRUE]
# colnames(Mido.data) <- read_tsv("~/scHAMR/STARsolo_results/Solo.out/Gene/filtered/barcodes.tsv", col_names=FALSE)[, 1, drop=TRUE]

# or:
# Mido.data <- ReadMtx(mtx ="~/scHAMR/STARsolo_results/Solo.out/Gene/filtered/matrix.mtx", cells="~/scHAMR/STARsolo_results/Solo.out/Gene/filtered/barcodes.tsv", features="~/scHAMR/STARsolo_results/Solo.out/Gene/filtered/features.tsv")

# setting up a Seurat object and displaying it
Mido <- CreateSeuratObject(counts = Mido.data, project = "Mido", min.cells = 3, min.features = 200)

Mido
```

2. Quality Control

```R
## Quality Control and cells selection

# percent of mitrochondrial genes
Mido[["percent.mt"]] <- PercentageFeatureSet(Mido, pattern = "^MT-")

# Violin plot for QC metrics
png("~/scHAMR/Seurat/qc_check/pre-qc_vlnplot.png", width = 800, height = 600, pointsize = 12)
VlnPlot(Mido, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# scatter plots for features
png("~/scHAMR/Seurat/qc_check/pre-qc_scatter.png", width = 800, height = 400, pointsize = 12)
plot1 <- FeatureScatter(Mido, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Mido, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# choosing cells with high quality (here, has more than 200 but less than 2500 reads and less than 5% mt genes)
Mido <- subset(Mido, subset= nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)

# rechecking Q
png("~/scHAMR/Seurat/qc_check/post-qc_vlnplot.png", width = 800, height = 600, pointsize = 12)
VlnPlot(Mido, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
		
png("~/scHAMR/Seurat/qc_check/post-qc_scatter.png", width=800, height=400, pointsize= 12)
plot1 <- FeatureScatter(Mido, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Mido, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()	
```

3. Preprocessing

```R
# Normalizing
Mido <- NormalizeData(Mido, normalization.method = "LogNormalize", scale.factor = 10000


# Features selection (select genes of high variability)
# 2000 genes of highest variability are selected
Mido <- FindVariableFeatures(Mido, selection.method = "vst", nfeatures = 2000)
# these are the top 10 of them to show in the plot below
top10 <- head(VariableFeatures(Mido), 10)
# plot the variable features
png("~/scHAMR/Seurat/features_selection/variable_features.png", width = 800, height = 400, pointsize = 12)
plot1 <- VariableFeaturePlot(Mido)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# Linear Transformation (Scaling)
all.genes <- rownames(Mido)
Mido <- ScaleData(Mido, features = all.genes)
```

4. Linear Dimensional Reduction with PCA

```R
# Calulcating PCA
Mido <- RunPCA(Mido, features = VariableFeatures(object = Mido))

# examining PCA results
print(Mido[["pca"]], dims = 1:5, nfeatures = 5) 

# visualizing PCA results in VizDimReduction(), DimPlot(), and DimHeatmap()
png("~/scHAMR/Seurat/PCA/VizDimLoadings.png", width = 800, height = 400, pointsize = 12)
VizDimLoadings(Mido, dims = 1:2, reduction = "pca")
dev.off()
png("~/scHAMR/Seurat/PCA/DimPlot.png", width = 800, height = 600, pointsize = 12)
DimPlot(Mido, reduction = "pca")
dev.off()
png("~/scHAMR/Seurat/PCA/DimHeatmap.png", width = 1200, height = 1200, pointsize = 12)
DimHeatmap(Mido, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

## Determine the dataset dimensionality
# determining the top PCs (and how many) that represent a robust compression for the dataset.
Mido <- JackStraw(Mido, num.replicate = 100, dims = 50)
Mido <- ScoreJackStraw(Mido, dims = 1:50)
# Visualizing the distribution of p-values for each PC with uniform distribution (the dashed line). Elbow plot can be used as an alternative.
png("~/scHAMR/Seurat/PCA/JackStrawPlot.png", width = 800, height = 400, pointsize = 12)
JackStrawPlot(Mido, dims = 1:50)
dev.off()
png("~/scHAMR/Seurat/PCA/ElbowPlot.png", width = 800, height = 600, pointsize = 12)
ElbowPlot(Mido)
dev.off()
```


### Option 2: Scanpy in Python

```python
# Python commands for clustering using Scanpy
import scanpy as sc
# ... additional commands for data processing and clustering
```

## 5. Running HAMR

1. Clusters

```bash

mkdir HAMR_clusters
cd clusters_BAM
# creating a for loop that runs over all BAM files in the directory clusters_BAM
for i in $(ls *)
do python2 ~/HAMRdirectory/HAMR-1.2/hamr.py $i  ~/scHAMR/reference_genome/<reference genome fasta file.fa>  models/euk_trna_mods.Rdata ~/scHAMR/HAMR_clusters  HAMR_$(basename $i)  30 10 0.05 H4 0.01 0.05 0.05
done
```

2. Cell-by-Cell (optional)

```bash

mkdir HAMR_cells
cd splitted_bams

# creating a for loop that runs over all BAM files in the directory clusters_BAM
for i in $(ls *)
do python2 ~/HAMRdirectory/HAMR-1.2/hamr.py $i  ~/scHAMR/reference_genome/<reference genome fasta file.fa>  models/euk_trna_mods.Rdata ~/scHAMR/HAMR_cells  HAMR_$(basename $i .pdf)  30 10 0.05 H4 0.01 0.05 0.05
done
```

3. Bulk (optional)

```bash

mkdir HAMR_Bulk
python2 ~/HAMRdirectory/HAMR-1.2/hamr.py ~/scHAMR/STARsolo_results/Aligned.sortedByCoord.out.bam  ~/scHAMR/reference_genome/<reference genome fasta file.fa>  models/euk_trna_mods.Rdata ~/scHAMR/HAMR_Bulk HAMR_results 30 10 0.05 H4 0.01 0.05 0.05

# for filtered reads according to the true cells detected by STARsolo
python2 ~/HAMRdirectory/HAMR-1.2/hamr.py ~/scHAMR/filtered_bam/filtered_bam  ~/scHAMR/reference_genome/<reference genome fasta file.fa>  models/euk_trna_mods.Rdata ~/scHAMR/HAMR_Bulk HAMR_results 30 10 0.05 H4 0.01 0.05 0.05
```

## 6. Example Runs

- Instructions and commands for example analyses on specific cell types.


## 7. Installing Prerequisites

1. STAR

Check the official resources here <https://github.com/alexdobin/STAR>

```bash
# get the latest STAR release
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xzf 2.7.11b.tar.gz
cd STAR-2.7.11b/source
make STAR

# or 
sudo apt install rna-star


# to run in the command line
STAR [options]
```

2. bamtools

```bash
sudo apt install bamtools

# to run in the command line
bamtools <commands>
```

3. samtools

```bash
sudo apt install samtools

# to run in the command line
samtools <commands>
```

4. SRA-toolkit

```bash
sudo apt install sra-toolkit
```

5. 10X Genomics subset-bam

This software is not officially supported by 10X Genomics. If you have troubles installing the software, run the following commands to install:

```bash
mkdir subset_bam
cd subset_bam 
wget https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux
mv subset-bam_linux subset-bam
chmod +x ./subset-bam
export PATH=$PATH:~/<full/path/to/subset-bam/directory>/subset_bam/
```
To operate the software, run one of the following commands:

```bash
./subset-bam  # if you are in the directory
subset-bam # if you exported the path
./ <full/path/to/subset-bam/directory>/subset-bam   # else
```

6. 10X Genomics Cell Ranger

Check the official resources here https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in

```bash
# download tar.gz file from the website or use your specific link you got after logging in above with the following command
# replace x.y.z with the specific version 
wget -O cellranger-<x.y.z>.tar.gz "<replace with your specific link>"

# unpack the downloaded file
tar -xzvf cellranger-x.y.z.tar.gz

# export the path
export PATH=~/<full path to the file location>/cellranger-x.y.z:$PATH

# to run in command line
cellranger <command>
```

7. R

```bash
# to install
sudo apt install r-base

# to check the version installed
R --version

# to run in command line
R
```

8. Python

```bash
# installing python 3
sudo apt install python3
sudo apt install python3-pip

# installing python 2
sudo apt install python2
```

9. Seurat Package in R

Check the official resources here https://satijalab.org/seurat/

```bash
# install system dependencies
sudo apt-get update
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
# run R in terminal
R
# install R dependencies
install.packages(c("curl", "openssl", "igraph", "httr", "leiden", "plotly"))
# install Seurat
install.packages("Seurat")
# close R 
q()
```

10. Scanpy Package in Python
Check the official resources here https://scanpy.readthedocs.io/en/stable/
```bash
# install scanpy
python3 -m pip install scanpy

# run python3 in terminal
python3
# check scanpy
import scanpy as sc
```
