#!/bin/bash  

module load bowtie2/
module load samtools/

## fasta files were taken from 
## https://www.ncbi.nlm.nih.gov/nuccore/M21017.1?report=fasta (drosophila)
## https://www.ncbi.nlm.nih.gov/nuccore/U13369.1?report=fasta (human)
## and saved in the below $DIR directory
DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/riboDNA"

CHIMERAGENOME="human-drosophila-rDNA"
CHIMERAGENOME_DIR="${DIR}/chimeras/${CHIMERAGENOME}" 
LOGS_DIR="${CHIMERAGENOME_DIR}/logs" 

if [ ! -d "$CHIMERAGENOME_DIR" ]; then
  mkdir -p  "$CHIMERAGENOME_DIR"
  mkdir -p $LOGS_DIR
fi

## combine fasta files 
cat ${DIR}/drosophila_rDNA_M21017.1.fa  ${DIR}/human_rDNA_U13369.1.fa > ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa

## create bowtie2 index
bowtie2-build  "${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa" "${CHIMERAGENOME_DIR}/${CHIMERAGENOME}" --threads 50 --verbose > "${LOGS_DIR}/bowtie2index-${CHIMERAGENOME}.log"

