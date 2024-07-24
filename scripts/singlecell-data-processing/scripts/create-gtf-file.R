
## Create GTF file

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Overview

## Use this script to create a GTF file for cellranger mkref
## Here I am using as example the chm13.draft_v2.0.gene_annotation.gff3 file from the T2T human genome and:
## 1. Convert it from gff3 to gtf
## 2. Extract only specific gene biotypes (ie, attributes)
#' 3. Remove duplicated genes (based on CAT and Liftoff gene annotation sources)
## Note that the original gff3 file was the UCSC GENCODEv35 CAT/Liftoff v2 obtained from the [T2T github webpage](https://github.com/marbl/CHM13?tab=readme-ov-file)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Set up the script

library(rtracklayer)
library(data.table)
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)

options(width=200)

gff.file = '/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/hs1/genes/chm13.draft_v2.0.gene_annotation.gff3'
gff <- import(gff.file)
gtf <- as.data.table(gff)

## Filter the GTF file to retain only those entries having gene biotypes of interest
kable(sort(table(copy(gtf)[,c('source_gene','gene_biotype')]%>%unique()%>%dplyr::pull('gene_biotype')),decreasing=T),align='c') ## this line prints in a visually pleasant format the number of entries in the gtf file per gene_biotype category

## select the biotypes you want based on the above output
biotypes = c(
    'protein_coding',
    'lncRNA',
    'IG_V_pseudogene',
    'IG_V_gene',
    'IG_D_gene',
    'IG_J_gene',
    'IG_J_pseudogene',
    'IG_C_gene',
    'IG_C_pseudogene',
    'IG_pseudogene',
    'TR_V_gene',
    'TR_V_pseudogene',
    'TR_J_gene',
    'TR_J_pseudogene',
    'TR_C_gene'
)

gtf.subset = copy(gtf)[gene_biotype %in% biotypes]

## In this GTF file, each gene name might have multiple sources of gene annotations (ie, CAT and Liftoff). 
## Gene IDs are instead unique and labelled as LOFF_* or CHM13_* depending on their source
## To remove redundancies (ie,many duplicated entries for a given gene name), if a gene name has multiple annotation sources keep only those with source = CAT; otherwise if a gene has only 1 annotation keep what's reported

gtf.subset = gtf.subset[,numberSources := length(unique(source)),by=.(source_gene_common_name)] ## this line counts how many annotation sources a given gene has
gtf.subset = gtf.subset[, keep := ifelse(numberSources==2, ifelse(source == 'CAT','keep','discard'),'keep')]
gtf.subset = gtf.subset[ keep=='keep']

## Print a table showing the number of genes belonging to each different annotation source (either CAT or Liftoff) in our GTF file
kable(sort(table(copy(gtf.subset)[,c('source_gene','source')]%>%unique()%>%dplyr::pull('source')),decreasing=T),align='c')

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Export files
gtf.subset <- makeGRangesFromDataFrame(gtf.subset,keep.extra.columns=T) ## reconvert data.frame to a GRange necessary for exporting

## location where you want to export the file. 
## Do not overwrite the file located within /stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/hs1/genes/ unless you have used the exact same above defined biotypes
output_gtf = '/vast/scratch/users/vespasiani.d/tmp/chm13.draft_v2.0.gene_annotation-gff3Converted-AttributesAndSourceFiltered.gtf' 
export(gtf.subset,output_gtf,"gtf") ## export info as a new gtf file...this might take a while

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()
