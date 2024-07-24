#!/bin/bash

## Set up ------------------------------------------------------------------

module load bcl2fastq/2.19.1 
module load pigz/
module load fastqc

## I/O variables to specify ------------------------------------------------------------------

RECEIVED_TAR_BCL="</path/to/tar/bcl/file>" ## this is the file you received straight from your sequencing provider
SAMPLESHEET="<path/to/sample/sheet/csv/file>" ## this is a csv file containing the lane, sampleID, and sequences of the single/dual indices used (ie, i7 and i5)
OUTDIR_BCL="</path/to/directory/where/to/extract/the/content/tarred/bcl>" ## this should be on vast/scratch 
OUTDIR_FASTQ="</path/to/directory/where/to/save/the/converted/fastq/files>" ## this is a directory that will contain all converted fastq files
OUTDIR_FASTQC="</path/to/where/you/want/to/store/the/fastQC/results>" 

## Extract bcl files  ------------------------------------------------------------------

tar -xvf "${RECEIVED_TAR_BCL}" --directory "${OUTDIR_BCL}"

## Convert bcl to fastq  ------------------------------------------------------------------

INTEROP_DIR="$OUTDIR_BCL/InterOp" ## leave this as is. It points to a sub-directory within the untarred bcl directory that bcl2fastq needs to be specified (perhaps check for its presence)

bcl2fastq  \
-R "${OUTDIR_BCL}"  \
--output-dir "${OUTDIR_FASTQ}" \
--interop-dir "${INTEROP_DIR}" \
--sample-sheet  "${SAMPLESHEET}" \
-r 10 -p 30 -w 10 \
--barcode-mismatches 1 \
--create-fastq-for-index-reads \
--minimum-trimmed-read-length 0 \
--mask-short-adapter-reads 0 \
--ignore-missing-positions \
--ignore-missing-controls \
--ignore-missing-filter \
--ignore-missing-bcls \
--use-bases-mask Y150n,I10,I10,Y150n ## change this to the your chosen read-sequencing configuration. This is for 150PE. Inspect the RunInfo.xml file within the untarred bcl dir to define that

# Run FastQC  ------------------------------------------------------------------

FASTQFILES="$(find $OUTDIR_FASTQ -type f -name \*.fastq.gz )" ## leave it as is. it will list all files ending with -name and run fastQC on them

## Create OUTDIR_FASTQC if not existing (required by fastQC)
if [ ! -d "$OUTDIR_FASTQC" ]; then
  mkdir -p  "$OUTDIR_FASTQC"
fi

fastqc --outdir ${OUTDIR_FASTQC} \
       --threads 40 \
       --nogroup \
       $FASTQFILES
