#!/usr/bin/env Rscript

## Author: davide vespasiani
## Mail: vespasiani.d@wehi.edu.au
## Description: Script to convert a gff3 file into a gtf 
## Script tested on R 4.4.0 

args = commandArgs(trailingOnly=TRUE)

packages <- rownames(installed.packages())

if (! 'rtracklayer' %in% packages) {
  BiocManager::install("rtracklayer")
}

if (length(args)==0) {
  stop("ERROR: you need to specify the full path + name of the gff3 file you want to convert", call.=FALSE)
}

gff_file = args[[1]]
gtf_output = paste(dirname(gff_file),paste(gsub('\\.gff3.*','',basename(gff_file)),'-gff3Converted.gtf',sep=''),sep='/')

cat('\n Importing the', basename(gff_file),' as a GFF file \n')

gff <- rtracklayer::import(gff_file)

cat('\n Exporting the converted', basename(gtf_output),'as a GFT file in the same input directory \n')

rtracklayer::export(gff,gtf_output,"gtf")

cat('\n Done with the GFF3-to-GTF conversion \n')


