# Calibrated PE Alignment
--------------------------------------------------

## Table of Contents
- [Calibrated PE Alignment](#calibrated-pe-alignment)
  - [Table of Contents](#table-of-contents)
  - [Project overview](#project-overview)
  - [Content](#content)
  - [Prerequisites](#prerequisites)
    - [Software dependencies](#software-dependencies)
  - [Quick start](#quick-start)
    - [Creating a metadata file](#creating-a-metadata-file)
    - [Running the pipeline on HPC](#running-the-pipeline-on-hpc)
  - [Output](#output)

--------------------------------------------------

## Project overview 

This repository contains all the necessary scripts to perform a calibrated whole-genome alignment of paired-end (PE) genomic bulk-sequencing data. Calibration here refers to the situation where you have added molecular spike-in DNA into your sample(s) before library preparation and you want to leverage them in order to improve the quantification and comparison of your actual sample(s). Assuming those spiked-in DNA molecules are not artificial but rather come from an organism (eg., E.coli) from which we have genome sequencing data already available, the idea is to create a chimera genome containing both the DNA sequences of your genome of interest (eg. hg38, mm19) and the spike-in genome (eg, ecoliASM584v2) and align your sequencing reads to that. <br/>

Specifically, the scripts contained within this repository will allow you to:

1. Create a chimera genome given 2 fasta files (1 for the spike-in genome and 1 for the genome of interest)
2. Align your PE sequencing reads to it 
3. Retain uniquely mapped (and optionally deduplicated) reads 
4. Extract the corresponding DNA fragments from the PE reads, either uniquely mapping to the spike-in genome or the genome of interest
5. Obtain a scaled genome coverage for the fragments mapping to the genome of interest

**NB:** All the scripts within this repository have been designed to be executed on the HPC. They have not yet been tested on a personal computer.

## Content

Below an overview of the content of this repository.

```
calibrated-pe-alignment/
├── config
│   └── modules.txt 
├── README.md
├── scripts
│   ├── calibrated-pe-alignment.sh
│   ├── create-chimera-genome.sh
│   └── parse-bowtie2-alignment-summary.R
|   └── combine-summary-tables.R
└── slurm
    ├── slurm-alignToChimera.sh
    └── slurm-createChimera.sh
```

* the `scripts/` folder contains all scripts necessary to execute the pipeline
* the `config/` folder contains a `.txt` file listing all modules (with versions) required to run the pipeline
* the `slurm/` folder contains two `.sh` files users can use to submit the pipeline either as an array or as a single sbatch command using SLURM on the HPC

--------------------------------------------------
## Prerequisites

To run this pipeline you fundamentally need to have 4 files:
1. Fasta file for the spike-in genome
2. Fasta file for the genome of interest
3. A sequencing fastq file(s) (one per read)
4. A 3-columns tab-separated file (check the [Creating a metadata file](#creating-a-metadata-file) section to see how to create it) containing:
   1. Sample name
   2. Name of sequencing read 1
   3. Name of sequencing read 2

### Software dependencies

To run the scripts within this repository you need to have installed:
1.[dplyr](https://dplyr.tidyverse.org/) R package in the same R module you will be using for running this pipeline. Here I am using R/4.4.0. So for example:

```
module load R/4.4.0 # load the module
R # enter an R session
install.packages('dplyr') # install the package. This might ask you from which CRAN repository you want to install it. Just select that corresponding to your geographical location
library(dplyr) # to ensure it has successfully installed
```

2.[pandas](https://pandas.pydata.org/), [sys](https://docs.python.org/3/library/sys.html), [openpyxl](https://openpyxl.readthedocs.io/en/stable/), [os](https://docs.python.org/3/library/os.html) python modules in the same python module you will be using for running this pipeline (same as above). Here I am using python/3.11.4. So for example:

```
module load python/3.11.4 # load the module
python ## enter a python session
import pandas,sys,os,openpyxl ## test is you already have these modules installed (some of them should).
```

If you dont get anything then it means you already have the modules installed in your python module. So you are good to go. Otherwise, if you instead get an error like:

```
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ModuleNotFoundError: No module named <moduleName>
```

Then it means you need to install that module. To do so you need to exit the python session (ie, `exit()`). Then from the command line run: `pip install <the-module-you-need> ## eg  pip install openpyxl`.

--------------------------------------------------

## Quick start

Below there is a step-by-step guide on how to run this pipeline. Before you run each script ensure you have your own space on the vast scratch file system on the Milton HPC. To check that run `cd /vast/scratch/users/$USER/` to go into that folder. 

**If you are a member of the Hamish King lab @ WEHI and if your chimera genome has already been created then go to 4.; otherwise go to 3.**

**If you already have a chimera genome directory containing a bowtie2 index and a 3 chrom.sizes.txt files (one per each individual genome and for the chimera) then skip to 4.**

1. Copy this repository into your project folder 
2. Ensure all scripts within the `calibrated-pe-alignment/scripts/` folder are executable and if not run `chmod +x /path/to/your/project/folder/calibrated-pe-alignment/scripts/*` from the command line
3. Create a chimera genome assembly by running:

```
/path/to/your/project/folder/calibrated-pe-alignment/scripts/create-chimera-genome.sh -s <spike-in-genome> -g <genome-of-interest>  -o <outputDir> 

## Flags description
# -s = full path to the fasta file (including filename) of the spike-in genome. File can be gzip compressed and it must end with .fa.gz
# -g = full path to the fasta file (including filename) of the genome of interest. File can be gzip compressed and it must end with .fa.gz
# -o = the path to the directory where you want to save all the output of this script. Here the script will create 2 directories: 1) named <spike-in-genome>-<genome-of-interest> containing all output files and 2) a logs directory containing the log from bowtie2-build
```
   
4. Align reads to the chimera genome by running:

```
/path/to/your/project/folder/calibrated-pe-alignment/scripts/calibrated-pe-alignment.sh -o <outDir> -d <dataDir> -g <chimeraGenomeDir>  -k <keepDuplicates> -m <metadataFile>

## Flags description
# -o = full path to the out directory where all final files and subdirectories will be saved
# -m = full path to the tab-separated metadata.txt file containing the sample information. The first 3 fields of this file should be: 1) sample name; 2) file name for forward read and 3) file name for reverse read. File must not have header
# -d = full path to the data directory containing your fastq files. Fastq files for the paired reads should be named the same as field 2 and field 3 of the metadata file
# -g = full path to the directory containing all files (ie, bowtie2 index, chromosome.sizes.txt) for the chimera genome previously generated by the create-chimera-genome.sh script. NB: file naming for the chimera genome should be <spikeGenome-genomeOfInterest>, eg ecoliASM584v2-hg38
# -k = yes/no depending on whether you want or not to keep duplicated reads
```

(Optional)
5. Combine the summary tables for each sample into a single file
  
```
Rscript /path/to/your/project/folder/calibrated-pe-alignment/scripts/combine-summary-tables.R /path/to/your/specified/output/directory/ filePattern

# For example
Rscript ./combine-summary-tables.R /your/specified/output/directory/ -summary.txt

```

This latter script will look for the `tables/` subdirectory within `/your/specified/output/directory/` and will take all samples inside that directory ending with your specified pattern (eg, *-summary.txt*). 
It will combine all files into a single table which will then be written at `/your/specified/output/directory/tables/all-samples-summary-table.txt`

### Creating a metadata file

The metadata file is required to run the `calibrated-pe-alignment.sh` script (ie to align fastq files to the chimera genome). This is a 3-columns no-header tab-separated txt file (see example below) containing the sample name (column1); name of fastq file containing the forward reads (column2) and name of fastq file containing the reverse reads (column3).  

```
S000488 G000454_batch1_1A       2024-03-20
S000488 G000454_batch1_1B       2024-03-20
S000488 G000454_batch1_2A       2024-03-20
S000488 G000454_batch1_2B       2024-03-20
```

There are many ways to create this file, either manually or programmatically. However, I noticed that uncorrect formatting of this file often causes unexpected failures of the main script. To avoid this, I've created another script (`scripts/prepare-metadata-file.py`) that takes as input a CSV or a xlsx file and generates a correct metadata file within the same directory of the input file. This script takes 4 arguments. Below you can see an example on how to run itt: 

```
path/to/your/project/folder/calibrated-pe-alignment/scripts/prepare-metadata-file.py  \
   csv/xlsx file \ ## this is the full path to the csv/xlsx file  
   R1 column ID \ ## the id of the column containing the fastq filenames for read 1 (ie forward read)
   R2 column ID  ## the id of the column containing the fastq filenames for read 2 (ie reverse read)
```
**Note that if using an xlsx file as input, the samples information you want to extract need to be on the first sheet.**

### Running the pipeline on HPC

To make it even easier I've provided 2 wrappers `.sh` scripts within the `slurm` directory that you can use to submit jobs on a SLURM cluster and:
1. create a chimera genome --> `slurm-createChimera.sh`
2. align fastq files to this genome --> `slurm-alignToChimera.sh`

To run them you just need to change:
* the `#SBATCH` configuration parameters
* the value of each variable within each script to ensure they correctly point to the location of your input files and to your output directory.

Note that `slurm-alignToChimera.sh` can be run either as an array or as a single job. I highly suggest you the former as it will submit *n* jobs with *n = the number of lines within the metadata file (which I remind you must not have any header)*. To run the job as an array you'll simply have to do: `sbatch --array=1-n slurm/slurm-alignToChimera.sh ## n = must be an integer`. Alternatively, you can always opt to run everything in a single instance by executing: `sbatch slurm/slurm-alignToChimera.sh`. <br/>

What will happen is that the `calibrated-pe-alignment.sh` script (wrapped within the provided slurm script) will loop through each line of the metadata file and produce output for each sample one after the other. The end results will be the same, it's just the timing to completion that will be quicker.

--------------------------------------------------

## Output

This pipeline will output all the final files into the user-specified output directory. However, all intermediate files will instead be written on the vast scratch filesystem on Milton in a newly created directory named as the main output directory. Both output directories (regardless whether containing intermediate or final files) will have the same structure but they'll differ in their file content. Below a outline of the content of the main output directory

```
main-output-directory
├── bam
│   ├── <your-sample>-fragments-uniqMap-<spikeInGenome>.bam  ## fragments uniquely mapped to the spike-in genome
│   ├── <your-sample>-fragments-uniqMap-<spikeInGenome>.bam.bai
│   ├── <your-sample>-fragments-uniqMap-<genomeOfInterest>.bam ## fragments uniquely mapped to the genome of interest
│   ├── <your-sample>-fragments-uniqMap-<genomeOfInterest>.bam.bai
├── bed
│   ├── <your-sample>-fragments-calibratedGenCov-<genomeOfInterest>.bg ## a bedGraph file containing a scaled genome coverage of the number of fragments uniquely mapping to the genome of interest (ie, scaled = calibrated on the number of fragments mapping to the spike-in genome) 
├── bigwig
│   ├── <your-sample>-fragments-calibratedGenCov-<genomeOfInterest>.bw  ## same as the bedGraph above but in a bigWig format for visualisation
├── logs
│   ├── <your-sample>-bowtie2-alignment-summary.txt ## a tabularised version of the log file generated by bowtie2 reporting the alignment stats
└── tables
    ├── <your-sample>-summary.txt # a table containing the alignment rates, percent duplicated reads (if decided to remove them), scaling factor, normalised library size and other stats for your sample
```

Now below instead you have an overview of the content of the temporary output directory

```
intermediate-output-directory
├── bam
│   ├── <your-sample>-fragments-uniqMap.bam  ## fragments uniquely mapped to the chimera genome
│   ├── <your-sample>-fragments-uniqMap-<spikeInGenome>.bam ## same file as the one on the main output dir but unsorted
│   ├── <your-sample>-fragments-uniqMap-<genomeOfInterest>.bam ## same as above
|   ├── <your-sample>-reads-uniqMap.bam  ## reads uniquely mapped to the chimera genome
│   ├── <your-sample>-reads-uniqMap-sortName.bam 
│   ├── <your-sample>-reads-uniqMap-sortName-fixmate.bam
│   ├── <your-sample>-reads-uniqMap-sortPos-fixmate.bam
│   ├── <your-sample>-reads-uniqMap-sortPos-fixmate-dedup.bam 
|   |── <your-sample>-reads-uniqMap-sortName-fixmate-final.bam ## a sorted bam by name with the same content of *sortPos-fixmate.bam or *sortPos-fixmate-dedup.bam depending on whether selected to keep or remove duplicated reads
├── bed
│   ├── <your-sample>-fragments-uniqMap.bed ## a be file containing a fragments uniquely mapping to the chimera genome
│   ├── <your-sample>-reads-uniqMap-sortName-fixmate-final.bed ## 10 column-bed file containing the paired-end reads delimiting the fragments
│   ├── <your-sample>-fragments-uniqMap-<genomeOfInterest>.bed ## same content as the above bam file but in a bed format 
│   ├── <your-sample>-fragments-uniqMap-<spikeInGenome>.bed ## same content as the above bam file but in a bed format 
├── bigwig
│   ├── <your-sample>-fragments-calibratedGenCov-<genomeOfInterest>.bw  ## same as the bedGraph above but in a bigWig format for visualisation
├── logs
│   ├── <your-sample>-bowtie2-alignment-summary.txt ## the original log file generated by bowtie2
└── tables
    ├── <your-sample>-summary.txt # an intermediate table containing the read/fragments counts for each file  
```

