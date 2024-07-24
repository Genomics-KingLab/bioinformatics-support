#!/bin/bash

## Set up ------------------------------------------------------------------

module load cellranger/8.0.0

## I/O variables to specify ------------------------------------------------------------------

GENOME="hs1" ## change this to whichever genome assembly you are using
GENOMES_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/${GENOME}"## path to the folder containing all info about your genome assembly
GTF_FILE="${GENOMES_DIR}/genes/chm13.draft_v2.0.gene_annotation-gff3Converted-AttributesAndSourceFiltered.gtf" ## leave it or change it to whichever gtf file you want to use to create your reference genome. Check readme for explanation on this file
OUTDIR="</path/where/to/save/cellranger/mkref/results>"

## Running cellranger mkref   ------------------------------------------------------------------

cellranger mkref \
 --genome="$GENOME" \
 --fasta="${GENOMES_DIR}/fasta/${GENOME}.fa" \
 --genes="${GTF_FILE}" \
 --nthreads=20 --memgb=100 \
 --output-dir "${OUTDIR}"
