## Single-cell data processing
--------------------------------------------------

## Table of Contents
- [Single-cell data processing](#single-cell-data-processing)
- [Table of Contents](#table-of-contents)
- [Project overview](#project-overview)
- [Content](#content)

--------------------------------------------------

## Project overview 

This folder contains a series of scripts to perform the main single-cell data processing steps. These include:
1. Bcl2fastq file conversion
2. Creating a reference genome for cellRanger
3. Aligning fastq files to the custom reference genome using cellRanger multi
4. Droplet processing which include sample demultiplexing and removal of low quality cells


## Content

Below I am listing the content of each subfolder contained within this repository 
1. `scripts/` --> all scripts necessary to process single-cell data
2. `examples/` --> a series of example files showing how information should be saved in order to run software (mostly cellranger)