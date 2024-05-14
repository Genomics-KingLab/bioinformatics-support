#!/usr/bin/env Rscript

## Author: Hamish King and Davide Vespasiani
## email: king.h@wehi.edu.au and vespasiani.d@wehi.edu.au
## Description: This script tabularise the log file produced by bowtie2 which contains all the alignment summary statistics
## Usage: Call this script from the command line specifying as arguments (in order):
## 1) the full path to the file containing the bowtie2 alignment summary statistics
## 2) the full path to the output directory where you want the table to be saved. The output file will be a tab-separated txt file named <your-sample>-bowtie2-alignment-summary.txt

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("ERROR: you need to specify the location of the bowtie2 log file and the path to the output directory as arguments", call.=FALSE)
}

if(!require("dplyr")){
  r <- getOption("repos")
  r["CRAN"] <- "http://cran.us.r-project.org"
  options(repos = r)
  install.packages("dplyr", dependencies =TRUE)
}

bowtie2Logfile = args[[1]]
outdir = args[[2]]
sample = gsub('.txt','',basename(bowtie2Logfile))
sample = gsub('-bowtie2-alignment-results','',sample)

alignRes = read.table(bowtie2Logfile,header=F,fill=T)

alignmentResults = data.frame(
        TotalPairedReads = as.numeric(as.character(alignRes$V1[2])),
        unMappedReads = as.numeric(as.character(alignRes$V1[3])),
        uniqMappedReads = as.numeric(as.character(alignRes$V1[4])),
        multiMappedReads = as.numeric(as.character(alignRes$V1[5])),
        overallReadsAlignmentRate = as.numeric(substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1))
)

alignmentResults <- mutate(alignmentResults, uniqueAlignmentRate = round(100 * (uniqMappedReads/TotalPairedReads),2))
alignmentResults <- mutate(alignmentResults, sample = sample,)
alignmentResults <- dplyr::select(alignmentResults,c('sample'),everything())

write.table(x = alignmentResults,file = paste(outdir, paste(sample,'-bowtie2-alignment-summary.txt',sep=''),sep='/'),quote=F,row.names=F,col.names=T,sep='\t')
