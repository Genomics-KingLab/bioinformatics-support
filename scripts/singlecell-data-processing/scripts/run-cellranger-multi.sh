#!/bin/bash

# Set up script ---------------------------------------------------

module load cellranger/8.0.0

## I/O variables to specify ------------------------------------------------------------------

SAMPLE="<name of the sample to analyse>"
OUTDIR="</path/where/to/save/all/cellranger/multi/output>"
CONFIG_FILE="</path/to/where/config/csv/file/is/located>" ## see examples/cellranger-multi-config.csv to see an example of how this file is formatted

## The lines below remove the output directory in case it already exists. This is required for cellranger multi to work
if [ -d "$OUTDIR" ]; then
  rm -rf $OUTDIR
fi

## Run cellranger multi  -----------------------------------------------------------------------------------------------------------------------------------

cellranger multi --id=${SAMPLE} \
                 --csv=${CONFIG_FILE} \
                 --localcores=20 \
                 --localmem=100 \
                 --output-dir=${OUTDIR} 
