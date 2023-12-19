# PhenoDriverR

R package for identifying personalized driver genes and exploring their oncogenic mechanism

## Installation
PhenoDriverR depends on org.Hs.eg.db, DESeq2, AnnotationDbi, clusterProfiler R packages from Bioconductor and please make sure they are installed before installing PhenoDriverR. An R version >= 4.1 is suggested.

To speed up matrix product, it is highly recommended to install openblas library first. For ubuntu, run the following command:
```shall
sudo apt install libopenblas-dev
```

Install the required packages
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "DESeq2", "AnnotationDbi", "clusterProfiler")
```

Install PhenoDriver
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("NWPU-903PR/PhenoDriverR")
```

## Data required
PhenoDriverR needs personalized gene expression and mutation data to identify driver genes. In addition, a weighted and directed Signal Tranduction Network (STN) also required for PhenoDriver workflow. Breast cancer data, as an example, can be found in [another github repository](https://github.com/NWPU-903PR/PhenoDriver-Paper)

## Usage of PhenoDriverR
See [the R Script](https://github.com/NWPU-903PR/PhenoDriver-Paper/blob/master/main.R) for an example of PhenoDriver full process execution.
