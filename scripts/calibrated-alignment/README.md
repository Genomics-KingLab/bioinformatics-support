# Calibrated PE Alignment
--------------------------------------------------

## Table of Contents
- [Calibrated PE Alignment](#calibrated-pe-alignment)
  - [Table of Contents](#table-of-contents)
  - [Project overview](#project-overview)
  - [Content](#content)
  - [Prerequisites](#prerequisites)
  - [Quick start](#quick-start)
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

To run this pipeline you need to fundamentally need to have 4 files:
1. Fasta file for the spike-in genome
2. Fasta file for the genome of interest
3. A sequencing fastq file(s) (one per read)
4. A metadata txt file. This is a 3-columns tab-separated file (with no header) containing in order:
   1. Sample name
   2. Name of sequencing read 1
   3. Name of sequencing read 2

You also need to have installed the [dplyr R package](https://dplyr.tidyverse.org/) in the same R module you will be using for running this pipeline. If this is R/4.4.0 then simply:

```
module load R/4.4.0 # load the module
R # enter an R session
install.packages('dplyr') # install the package. This might ask you from which CRAN repository you want to install it. Just select that corresponding to your geographical location
library(dplyr) # to ensure it has successfully installed
```

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
│   ├── <your-sample>-fragments-uniqMap-<genomeOfInterest>.bed ## same content as the above bam file but in a bed format --> to be moved to tmp folder
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
├── bigwig
│   ├── <your-sample>-fragments-calibratedGenCov-<genomeOfInterest>.bw  ## same as the bedGraph above but in a bigWig format for visualisation
├── logs
│   ├── <your-sample>-bowtie2-alignment-summary.txt ## the original log file generated by bowtie2
└── tables
    ├── <your-sample>-summary.txt # an intermediate table containing the read/fragments counts for each file  
```

