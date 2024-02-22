# scHAMR Pipeline Manual

## Table of Contents

1. [scHAMR pipeline diagram](#1-pipeline-diagram)
2. [Prerequisites](#2-prerequisites)
3. [Preprocessing](#3-preprocessing)
   - [Starting up and loading the dataset](#31-starting-up-and-loading-the-dataset)
   - [Building the genome index](#32-building-the-genome-index)
   - [Alignment and Quantification](#33-alignment-and-quantification)
4. [Processing](#4-processing)
   - [Cell by cell](#41-cell-by-cell-analysis-optional)
   - [Clustering](#42-clustering)
     - [Option 1: Seurat in R](#option-1-seurat-in-r)
     - [Option 2: Scanpy in Python](#option-2-scanpy-in-python)
   - [Generating cluster-specific BAM files](#43-generating-cluster-specific-bam-files)
5. [Running HAMR](#5-running-hamr)
6. [Example Runs](#6-example-runs)
   - [Drosophila Escort Cells](#drosophila-escort-cells)
   - [Human Pancreatic Islets](#human-pancreatic-islets)
7. [Installing Prerequisites](#7-installing-prerequisites)
8. [Troubleshooting]()

## 1. Pipeline Diagram
<details>
<img src="./pip1.png" alt="alt text" width="500" height="300">

<img src="./pip2.png" alt="alt text" width="500" height="300">

</details>

## 2. Prerequisites
<details>
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
</details>

## 3. Preprocessing
<details>

### 3.1. Starting up and Loading the Dataset

<details>
Commands to set up directories, load, and preprocess data

```bash
# creating a directory that will hold all the analyses
mkdir -p scHAMR
cd ~/scHAMR

# loading the sample data from GEO in SRA format
mkdir -p SRR_data
cd SRR_data
prefetch <SRRxxxxxxxx> --max-size 200G
cd ~/scHAMR

# converting it to FASTQ and spliting the files of reads (R1, R2, and possibly I1...).
mkdir -p FASTQ_data
cd FASTQ_data
fasterq-dump ~/scHAMR/SRR_data/<SRRxxxxxxxx> --split-files
ls
cd ~/scHAMR
```
</details>

### 3.2. Building the Genome Index
<details>
Although the annotations should not be added in building the genome index since HAMR requires no spliced junctions, STARsolo requires the annotations to run and produce the count matrix after aligning. The spliced junction problems will be solved during the aligning step. Additionally, the annotations file needs to be filtered for exons as recommended by 10Xgenomics and STARsolo to properly create the count matrix.


Commands for building genome index with STAR
 
```bash
# Ensuring the starting directory is ~/scHAMR
cd ~/scHAMR

# genome index directory
mkdir -p reference_genome
cd reference_genome

# loading required files: genome fasta and annotations GTF
wget <link for ensembl reference genome fasta file>
wget <link for ensembl annotations GTF file>
gzip -d *.gz

# filtering the annotations for exons
cellranger mkgtf <input.annotations_file.gtf> <output.annotations_filtered_file.gtf> --attribute=gene_biotype:protein_coding

# building genome index using STAR
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir STAR_annotated_index/ --genomeFastaFiles <reference genome file.fa> --sjdbGTFfile <annotations_filtered_file.gtf> --genomeSAindexNbases 12 --genomeSAsparseD 3
cd ~/scHAMR
```
</details>


### 3.3. Alignment and Quantification
<details>

STARsolo is used for this step since it provides flexibility in use to work within the constraints of HAMR as well as producing comparable results to CellRanger
1.	Spliced junctions for mRNA need to be filtered out. To filter them out in bulk RNA-seq, the STAR aligner parameter --alignIntronMax 1 is usually used along with not including the annotations file in the genome index step. That is because --alignIntronMax 1 only controls the unannotated junctions and has no control over the junctions annotated in building the genome index. However, the annotations file is required for scRNA-seq as discussed earlier in the genome index step. To fix this problem, the annotated spliced junctions can be filtered out by increasing increasing the overhang to a number bigger than the read length, --alignSJDBoverhangMin 999 (n> read length).
2.	The “CB” tag must be included in the --outSAMattributes and that the produced file is a sorted BAM (--outSAMtype BAM SortedByCoordinate ) because the "CB" tag will not be included otherwise and a sorted BAM is also a requirement by HAMR. The “CB” tag will be used to generate BAM file for each cell or cluster later.
3.	HAMR requires only uniquely mapped reads. The parameter --outFilterMultimapNmax 1 is used to filter out multiple mapped reads.
4.	Some parameters such as –soloType, --soloUMIlen and input fastq files are adjusted acording to the used kit. For example. here, the parameters are adjusted for the 10x chromium 3" V2 kit. For the 10X chromium 3" V3, add --soloUMIlen 12. Additionally, STARsolo requires the 10x Genomics cells barcodes whitelist, which is different for different kit versions, to check for correct CBs. Review the STAR Aligner manual for more details and guidance.
5.	HAMR and Subset-bam require the BAM to be sorted and indexed.
6. Find the 10X barcodes whitelist [here](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html)


Commands for aligning the Reads, CB Demultiplexing, UMI Deduplication, Counting and Cell Calling with STAR

```bash
# loading the 10x Genomics cells barcodes whitelist
mkdir -p CB_whitelist
cd CB_whitelist
wget <link for CB whitelist txt file>
gzip -d *.gz
cd ~/scHAMR

# mapping with STARsolo
STAR --runThreadN 4   --genomeDir reference_genome/STAR_annotated_index/ --readFilesIn FASTQ_data/<Second file with actual cDNA reads.fastq>   FASTQ_data/<first file with CB(16b)+UMI(10b) reads.fastq>  --outFileNamePrefix STARsolo_results/   --outReadsUnmapped Fastx   --outSAMattributes NH   HI   NM   MD  CB UB sM sS sQ    --outFilterMultimapNmax 1   --outFilterMatchNmin 30   --outFilterMismatchNmax 4   --alignIntronMax 1   --alignSJDBoverhangMin 999   --soloType CB_UMI_Simple --soloCellFilter EmptyDrops_CR  --soloCBwhitelist CB_whitelist/<CB whitelist file.txt> --soloBarcodeReadLength 1 --soloCBlen 16 -- soloUMIlen <10 or 12 based on the 10X version> --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000

# indexing the resulted BAM
samtools index ~/scHAMR/STARsolo_results/Aligned.sortedByCoord.out.bam
```

The output of STARsolo includes the BAM file as well as raw and filtered count matrix in addition to other complementary files as summaries and logs. The filtered count matrix and BAM are required for the next steps.
</details>

</details>

## 4. Processing

<details>

### 4.1. Cell by Cell Analysis (Optional)

<details>

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

</details>

### 4.2. Clustering

<details>

``` bash
# directory for all clustering analyses
mkdir clustering
cd clustering
```

### Option 1: Seurat in R

<details>

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

</details>

### Option 2: Scanpy in Python

<details>

**1. Setting up and Loading dataset**

Scanpy expects zipped files:
```bash
# gzip STARsolo results for Scanpy
mkdir ~/HAMR/STARsolo_results/Solo.out/Gene/filtered/gzipped
gzip -c ~/HAMR/STARsolo_results/Solo.out/Gene/filtered/* ~/HAMR/STARsolo_results/Solo.out/Gene/filtered/gzipped/

# start in the ~/scHAMR directory
cd ~/scHAMR

# starting the python environment
python3
```

Importing libraries:
```python
import scanpy as sc
import os
import anndata
import scipy as sp
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.sparse import csr_matrix
from sklearn.model_selection import StratifiedShuffleSplit, train_test_split
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
np.random.seed(223)
```

Loading dataset:
```python
# Reading the dataset into python as an anndata
files_path = '../scHAMR/STARsolo_results/Solo.out/Gene/filtered/gzipped/'
adata = sc.read_10x_mtx(files_path)

# Displaying the AnnData object description
adata
```

**2. Data Cleaning and Quality Control**

In case there are multiple samples in the dataset, the cleaning needs to be done to individual samples.

Subsetting the dataset for individual samples:
```python
# Subset the data for each of the 5 samples based on the sample annotations in the data
unique_samples = adata.obs['sample'].unique()
sample_data = {}
for sample in unique_samples:
    sample_data[sample] = adata[adata.obs['sample'] == sample].copy()
    
#adata_<sample1_ID> = sample_data['<sample1_ID>']
#adata_<sample2_ID> = sample_data['<sample2_ID>']
#adata_<sample3_ID> = sample_data['<sample3_ID>']
#adata_<...> = sample_data['<...>']
```

Define a function for quality control check using metrics and visualizations:
```python
def data_quality_control_check(adata):
    """
    Perform quality control (QC) analysis on an AnnData object used in scRNA-seq data analysis.

    This function calculates and adds QC metrics to the AnnData object for each cell. These metrics include 
    the total counts of RNA molecules per cell, the number of detected genes, and the fraction of 
    mitochondrial (MT) genes. It also generates violin plots and scatter plots for these metrics to assist 
    in determining appropriate threshold values for further quality control filtering.

    Parameters:
    -----------
    adata : AnnData
        An AnnData object containing scRNA-seq data. This object should have cells as rows and genes as columns.
    
    Returns:
    --------
    adata : AnnData
        The modified AnnData object with added QC metrics. The metrics added are 'n_genes_by_counts' (number of 
        genes detected in each cell), 'total_counts' (total number of RNA molecules detected per cell), and 
        'pct_counts_mt' (percentage of counts belonging to mitochondrial genes).
    """

    # Identify and annotate mitochondrial genes, which start with MT in their ID
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    # Calculate quality check metrics, particularly: total counts, no. of genes, and MT genes fraction
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
    
    # Produce a violin plot for the quality check metrics 
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
    
    # Produce scatter plots for total count vs mitochondrial genes and gene count
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    
    return adata
```

Define a function that applies selected quality control metrics:
```python
def data_quality_control_apply(adata, min_counts, max_counts, min_genes, max_genes, max_pct_mt):
    
    """
    Apply a series of quality control filters to an AnnData object from scRNA-seq data.

    This function performs several filtering steps to remove low-quality cells based on specified 
    criteria: the maximum total counts, the minimum and maximum number of genes expressed, and the 
    maximum percentage of mitochondrial gene counts. The function prints the number of cells in the 
    dataset after each filtering step for tracking the impact of each criterion.

    Parameters:
    -----------
    adata : AnnData
        An AnnData object containing single-cell RNA sequencing data, with cells as rows and genes as columns.
    max_counts : int
        Maximum allowed total counts (sum of all gene expression counts) per cell. Cells exceeding this 
        threshold will be filtered out.
    min_genes : int
        Minimum number of genes that must be expressed in a cell. Cells with fewer expressed genes will 
        be filtered out.
    max_genes : int
        Maximum number of genes that must be expressed in a cell. Cells with more expressed genes will 
        be filtered out.
    max_pct_mt : float
        Maximum allowed percentage of mitochondrial gene counts. Cells with a higher percentage will be 
        filtered out.

    Returns:
    --------
    AnnData
        The filtered AnnData object.
    """
    
    # Number of cells before any filtering
    print('Total number of cells before filtering: {:d}'.format(adata.n_obs))
    
    # Filter out counts over min_counts
    sc.pp.filter_cells(adata, min_counts = min_counts)
    print('Number of cells after min count filter: {:d}'.format(adata.n_obs))
    
    # Filter out counts over max_counts
    sc.pp.filter_cells(adata, max_counts = max_counts)
    print('Number of cells after max count filter: {:d}'.format(adata.n_obs))

    # Filter out cells with under min_genes genes
    sc.pp.filter_cells(adata, min_genes = min_genes)
    print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

    #Filter out cells with over max_genes genes
    sc.pp.filter_cells(adata, max_genes = max_genes)
    print('Number of cells after gene filter: {:d}'.format(adata.n_obs))
    
    # Filter out cells with high percentage of mitochondrial genes
    #adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :].copy()
    adata = adata[adata.obs['pct_counts_mt'] < max_pct_mt].copy()

    print('Number of cells after MT pct filter: {:d}'.format(adata.n_obs))
    
    return adata
```

Repeat the following as neccessary for all samples individually:
```python
# Calculate and visualize QC metrics
data_quality_control_check(adata)

# Apply selected QC values
adata = data_quality_control_apply(adata61, min_counts=200, max_counts=50000, min_genes=1000, max_genes=5000, max_pct_mt=20)
adata

# Calculate and visualize QC metrics after applying QC
data_quality_control_check(adata)
```

If applicable, integrate all the samples back into one dataset:
```python
# Concatenate all the samples in adata
adata = anndata.concat(
    {bc: ad for bc, ad in zip(['<sample1_ID>', '<sample2_ID>', '<sample3_ID>', '<...>'], [adata_<sample1_ID>, adata_<sample2_ID>, adata_<...>])},
    label='sample',
    merge="same"
)
adata
```

Filtering out genes expressed in less than X cells in the whole dataset. Those genes are probably artifacts:
```python
# Filter genes
print('Total number of genes: {:d}'.format(adata.n_vars))

# Basic filtering - a Gene must be at least detected in 30 cells
sc.pp.filter_genes(adata, min_cells=30)

print('Number of genes after filtering: {:d}'.format(adata.n_vars))
```

**3. Data Transformation**

Normalizing:
```python
# Keep the count data that is not normalized in a counts layer.
adata.layers["counts"] = adata.X.copy()

# Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell.
sc.pp.normalize_total(adata, target_sum=1e4)
```

Logarithmizing:
```python
# Logarithmize the data
sc.pp.log1p(adata)

# Save the normalized and logarithmized raw data in the .raw attribute of the anndata.
adata.raw = adata
```

**4. Extracting Highly Variable Genes (HVGs) and Further Cleaning**

Highly variable genes analysis:
```python
# Extracting highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Visualization of highly variable genes
sc.pl.highly_variable_genes(adata)

# Get only highly variable genes
adata = adata[:, adata.var.highly_variable]
```

Regress out the effects of confounding variables, such as total count and mitochondrial genes percentage:
```python
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
```

Scaling the data to unit variance and clip values above standard deviation of 10 to minimize the effects of the outliers
```python
sc.pp.scale(adata, max_value=10)
```

**5. Batch Correction**

- Batch effects are systematic differences in data that arise not from biological variations but from technical or experimental differences when having multiple samples in the dataset. Batch correction in scRNA-seq is a critical step for ensuring that subsequent analyses reflect true biological differences rather than technical artifacts. A commonly used method is Combating Batch Effects (ComBat), which uses an empirical Bayes framework to model and adjust for batch effects in gene expression data. This method leverages the strengths of both frequentist (data-driven estimation) and Bayesian (prior and posterior updating) methods, making it particularly effective for high-dimensional genomic data where direct parameter estimation might be challenging.

Apply batch correction if the dataset has multiple samples:
```python
# Batch Correction using ComBat
sc.pp.combat(adata, key='sample')
```

**6. Principal Component Analysis (PCA)**

```python
# Performing Principal Component Analysis (PCA)
sc.tl.pca(adata, n_comps=100, svd_solver='arpack')

# Visualize elbow plot for the variance ratio across principal components
sc.pl.pca_variance_ratio(adata, n_pcs=100, log=True)
```

**7. Neighborhood Graph Construction**

```python
# Calculate the neighborhood graph with 50 PCs. 
sc.pp.neighbors(adata_explore, n_pcs=50, n_neighbors = 15)
```

**8. Uniform Manifold Approximation and Projection (UMAP)**

```python
# UMAP calculation
sc.tl.umap(adata_explore) # min_dist = 0.5 by default

# Visualize the UMAP with highlighting the different five samples in the dataset
sc.pl.scatter(adata_explore, basis='umap', color='sample')
```

**9. Clustering with Louvain or Leiden and Resolution Tuning**

The Louvain model or Leiden model can be selected here. This is an example with Leiden. The process is exactly the same with Louvain. Since a better performance is usually dataset specific, comparing both models is recommended.

Running clustering with different resolution parameter values:
```python
# Define a range of resolution values
resolutions = np.arange(0.1, 2.1, 0.1)  

# Initialize dictionaries to store scores
silhouette_scores = {}
davies_bouldin_scores = {}
calinski_harabasz_scores = {}

for res in resolutions:
    # Perform Louvain clustering at the given resolution
    sc.tl.leiden(adata, resolution=res, key_added=f'clustering_{res}')

    # Retrieve the cluster labels
    labels = adata.obs[f'clustering_{res}']

    # Assuming adata.obsm['X_pca'] contains the PCA reduced data
    X_pca = adata.obsm['X_pca']

    # Calculate Silhouette Score
    silhouette_scores[res] = silhouette_score(X_pca, labels)

    # Calculate Davies-Bouldin Score
    davies_bouldin_scores[res] = davies_bouldin_score(X_pca, labels)

    # Calculate Calinski-Harabasz Score
    calinski_harabasz_scores[res] = calinski_harabasz_score(X_pca, labels)

    # Plot the clusters
    sc.pl.scatter(adata, basis='umap', color=f'clustering_{res}', title=f'Resolution {res}')
```

Inspect the statistical metrics and plot them against resolution parameters:
```python
# Inspect the scores in the dictionaries
print("Silhouette Scores:", silhouette_scores)
print("Davies-Bouldin Scores:", davies_bouldin_scores)
print("Calinski-Harabasz Scores:", calinski_harabasz_scores)

# Plotting Silhouette Score
plt.figure(figsize=(10, 6))
plt.plot(list(silhouette_scores.keys()), list(silhouette_scores.values()), marker='o')
plt.xlabel("Resolution")
plt.ylabel("Silhouette Score")
plt.title("Silhouette Score for Different Resolutions in Louvain Clustering")
plt.show()

# Plotting Davies Bouldin Score
plt.figure(figsize=(10, 6))
plt.plot(list(davies_bouldin_scores.keys()), list(davies_bouldin_scores.values()), marker='o')
plt.xlabel("Resolution")
plt.ylabel("Davies Bouldin Score")
plt.title("Davies Bouldin for Different Resolutions in Louvain Clustering")
plt.show()

# Plotting Calinski Harabasz Score
plt.figure(figsize=(10, 6))
plt.plot(list(calinski_harabasz_scores.keys()), list(calinski_harabasz_scores.values()), marker='o')
plt.xlabel("Resolution")
plt.ylabel("Calinski Harabasz Score")
plt.title("Calinski Harabasz for Different Resolutions in Louvain Clustering")
plt.show()
```
After deciding on a resolution parameter, visualize the clustering results:
```python
sc.tl.leiden(adata, key_added= 'clusters', resolution=0.4)
sc.pl.scatter(adata, basis='umap', color=['clusters'])

labels = adata.obs['clusters']
X_pca = adata.obsm['X_pca']

# Calculate Silhouette Score
silhouette_score_value_leiden = silhouette_score(X_pca, labels)
print("Silhouette score: ", silhouette_score_value_leiden)

# Calculate Davies-Bouldin Score
davies_bouldin_score_value_leiden = davies_bouldin_score(X_pca, labels)
print("Davies Bouldin score: ",  davies_bouldin_score_value_leiden)

# Calculate Calinski-Harabasz Score
calinski_harabasz_score_value_leiden = calinski_harabasz_score(X_pca, labels)
print("Calinski Harabasz score: ", calinski_harabasz_score_value_leiden)
```

**10. (Optional) Assessing Clustering Robustness on Dataset Subsets**

Perform stratified sampling:
```python
# Define bins or regions for splitting

# Get UMAP coordinates
umap_coords = adata.obsm['X_umap']

# Define grid boundaries (these could be based on quantiles or other criteria)
x_bins = np.linspace(min(umap_coords[:,0]), max(umap_coords[:,0]), 100)
y_bins = np.linspace(min(umap_coords[:,1]), max(umap_coords[:,1]), 100)

# Digitize the UMAP coordinates to bin indices
x_bin_indices = np.digitize(umap_coords[:,0], x_bins)
y_bin_indices = np.digitize(umap_coords[:,1], y_bins)

# Combine the bin indices to form a stratification key
stratification_key = x_bin_indices * (1000) + y_bin_indices
```

```python
# Splitting the data into two subsets. Here, we perform stratified splitting. 
# However, if a bin has only one sample, we randomly assign the sample to a subset.

# Check if any bin has fewer than 2 samples
unique, counts = np.unique(stratification_key, return_counts=True)
# Find bins with fewer than 2 samples
bins_with_fewer_than_two = unique[counts == 1]

# Create masks for cells in bins with at least two samples and with only one sample
mask_fewer_than_two = np.isin(stratification_key, bins_with_fewer_than_two)
mask_at_least_two = ~mask_fewer_than_two

# Split indices into two groups
indices_fewer_than_two = np.where(mask_fewer_than_two)[0]
indices_at_least_two = np.where(mask_at_least_two)[0]

# Perform stratified split on cells in bins with at least two samples
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.5, random_state=42)
for first_subset_idx_stratified, second_subset_idx_stratified in sss.split(X=np.zeros(
    len(indices_at_least_two)), y=stratification_key[mask_at_least_two]):
    pass

# Adjust indices to original data
first_subset_idx_stratified = indices_at_least_two[first_subset_idx_stratified]
second_subset_idx_stratified = indices_at_least_two[second_subset_idx_stratified]

# Randomly assign cells from bins with fewer than two samples
first_subset_idx_random, second_subset_idx_random = train_test_split(indices_fewer_than_two, test_size=0.5, random_state=42)

# Combine indices from both stratified and random splits
first_idx = np.concatenate((first_subset_idx_stratified, first_subset_idx_random))
second_idx = np.concatenate((second_subset_idx_stratified, second_subset_idx_random))

# Create First and Second subsets
first_subset_adata = adata[first_idx].copy()
second_subset_adata = adata[second_idx].copy()
```

```python
print(first_subset_adata,'\n \n', second_subset_adata)
```


Running clustering on the dataset subsets:
```python
# First subset:
sc.tl.leiden(first_subset_adata, key_added= 'clusters_fsubset', resolution=0.4)
sc.pl.scatter(first_subset_adata, basis='umap', color=['clusters_fsubset'])

labels_first_subset = first_subset_adata.obs['clusters_fsubset']
X_pca_first_subset = first_subset_adata.obsm['X_pca']

# Calculate Silhouette Score
first_subset_silhouette_score = silhouette_score(X_pca_first_subset, labels_first_subset)
print("First Subset's Silhouette score: ", first_subset_silhouette_score)

# Calculate Davies-Bouldin Score
first_subset_davies_bouldin_score = davies_bouldin_score(X_pca_first_subset, labels_first_subset)
print("First Subset's Davies Bouldin score: ",  first_subset_davies_bouldin_score)

# Calculate Calinski-Harabasz Score
first_subset_calinski_harabasz_score = calinski_harabasz_score(X_pca_first_subset, labels_first_subset)
print("First Subset's Calinski Harabasz score: ", first_subset_calinski_harabasz_score)
```

```python
# Second subset:
sc.tl.leiden(second_subset_adata, key_added= 'clusters_ssubset', resolution=0.4)
sc.pl.scatter(second_subset_adata, basis='umap', color=['clusters_ssubset'])

labels_second_subset = second_subset_adata.obs['clusters_ssubset']
X_pca_second_subset = second_subset_adata.obsm['X_pca']

# Calculate Silhouette Score
second_subset_silhouette_score = silhouette_score(X_pca_second_subset, labels_second_subset)
print("Second Subset's Silhouette score: ", second_subset_silhouette_score)

# Calculate Davies-Bouldin Score
second_subset_davies_bouldin_score = davies_bouldin_score(X_pca_second_subset, labels_second_subset)
print("Second Subset's Davies Bouldin score: ",  second_subset_davies_bouldin_score)

# Calculate Calinski-Harabasz Score
second_subset_calinski_harabasz_score = calinski_harabasz_score(X_pca_second_subset, labels_second_subset)
print("Second Subset's Calinski Harabasz score: ", second_subset_calinski_harabasz_score)
```

**11. Exporting a Dataframe for Cell Barcodes and Clusters ID**

```python
# Extracting cell barcodes and cluster IDs
cell_barcodes = adata.obs_names
cluster_ids = adata.obs['clusters']

# Creating a dataframe
cb_cluster_df = pd.DataFrame({'CellBarcode': cell_barcodes, 'ClusterID': cluster_ids})

# Saving the dataframe to CSV
csv_file_path = "~/scHAMR/clustering/CBs_Clusters_dataframe.csv"  # Adjust the path as needed
cb_cluster_df.to_csv(csv_file_path, index=False)
```

**12. (Optional) Gene Markers and Cell Typing Analysis and Visualization**

Identifying significantly enriched genes in each cluster:
```python
# t-test
sc.tl.rank_genes_groups(adata, 'louvain_clusters', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```

Visualizing the expression of specific marker genes known for cell types on the UMAP:
```python
sc.pl.scatter(adata, basis='umap', color=['Gene_Marker_ID1', 'Gene_Marker_ID2', '...'])
```

Dot plot visualization for gene markers expression in clusters and corresponding cell type:

```python
# Marker genes to cell type dictionary
marker_genes_dict = {
    'cell_type1': ['Gene_Marker_ID1'],
    'cell_type2': ['Gene_Marker_ID2'],
    <....>
    'cell_type3': ['Gene_Marker_ID3', 'Gene_Marker_ID4'],
}

# Dotplot for visualization
sc.pl.dotplot(adata, marker_genes_dict, 'clusters')
```

Cell-type annotations:
```python
# Using the previous information to annotate cell types

map_names = {}

for c in adata.obs['louvain_clusters'].cat.categories:
    if c in ['0','4', '1', '2']:
        map_names[c] = 'cell_type1'
    elif c in ['5', '3']:
        map_names[c] = 'cell_type2'  
    elif c in ['6']:
        map_names[c] = 'cell_type3'
    elif c in ['<....>']:
        map_names[c] = 'cell_type<...>'
    else:
        map_names[c] = c

adata.obs['clusters_annotations'] = adata.obs['clusters']
adata.obs['clusters_annotations'] = adata.obs['clusters_annotations'].map(map_names).astype('category')
adata.obs['Cell_Type'] = adata.obs['clusters_annotations'].cat.reorder_categories(
    ['cell_type1', 'cell_type2', 'cell_type3', 'cell_type<...>'])
```

Finally, exit python environment

```python
exit()
```

</details>

</details>

</details>

### 4.3. Generating cluster-specific BAM files


## 5. Running HAMR
<details>

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
</details>

## 6. Example Runs
<details>

- Instructions and commands for example analyses on specific cell types.
<details>

### Drosophila Escort Cells 

</details>

<details>

### Human Pancreatic Islets

</details>

</details>

## 7. Installing Prerequisites

<details>

1. STAR

Check the official documentation [here](https://github.com/alexdobin/STAR>)

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

This software is not officially supported by 10X Genomics. Check the documentation [here](https://github.com/10XGenomics/subset-bam)

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

Check the official documentation [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in)

```bash
# download tar.gz file from the website or use your specific link you got after logging in above with the following command
# replace x.y.z with the specific version 
wget -O cellranger-<x.y.z>.tar.gz "<replace with your specific link>"

# unpack the downloaded file
tar -xzvf cellranger-<x.y.z>.tar.gz

# export the path
export PATH=~/<full path to the file location>/cellranger-<x.y.z>:$PATH

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

Check the official documentation [here](https://satijalab.org/seurat/)

```bash
# install system dependencies
sudo apt-get update
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libglpk-dev
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
Check the official documentation [here](https://scanpy.readthedocs.io/en/stable/)
```bash
# install scanpy
python3 -m pip install scanpy

# run python3 in terminal
python3
# check scanpy
import scanpy as sc
```
</details>

## 8. Troubleshooting
<details>
</details>