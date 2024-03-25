# Running the Scripts

Ensure the working directory is the directory that contains these files (scripts) or modify the paths as necessary.

### Starting up and Loading the Data

```bash
chmod +x ./startingup-loadingdata.sh
./startingup-loadingdata.sh
```

### Reference Genome

```bash
chmod +x ./reference-genome.sh
./reference-genome.sh
```

### Alignment and Quantification

```bash
chmod +x ./alignment-quantification.sh
nohup ./alignment-quantification.sh &
```

### Clustering

Run the Jupyter notebook clustering.ipynb after:

```bash
mkdir ~/scHAMR/clustering
```

### Generating clusters-specific bam files
```bash
chmod +x ./clusters-bam.sh
nohup ./clusters-bam.sh &
```


### Running HAMR
```bash
chmod +x ./running-hamr.sh
nohup ./running-hamr.sh &
```

