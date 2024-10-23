## Content

This folder contains a list of simplified scripts to analyse scRNA-seq data from the raw BCL data using CellRanger and Seurat v5 workflows.

Preliminary steps:
* Create a cellranger reference genome for your assembly of interest
* Create a GTF file containing the gene information you want for the genome of interest

Analysis:
1. Convert sequencing data from bcl to fastq
2. Align fastq reads to reference genome using cellranger multi in case your fastq files contain a combination of Gene Expression, Feature Barcode, and V(D)J libraries generated from a single GEM well, otherwise use cellranger count
3. Droplet processing of the cellranger output which includes sample (eg, donor) demultiplexing ([HTODemux](https://satijalab.org/seurat/articles/hashing_vignette.html)) and creation of Seurat object containing the count matrices for the different assays
4. Single-cell data processing using Seurat v5 workflow (check [installation](https://satijalab.org/seurat/articles/install_v5) and [PBMC tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial))
