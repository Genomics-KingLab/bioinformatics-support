#!/usr/bin/env Rscript

## Author: Hamish King and Davide Vespasiani
## email: king.h@wehi.edu.au and vespasiani.d@wehi.edu.au
## Description: This script greps a list of tab-delimited files located within the inputDir and ending with the specified filePattern and combines them into a single summary table which will be exported into the same inputDir 
## Usage: Call this script from the command line specifying as arguments (in order):
## 1) inputDir --> the full path to the input directory containing all the tab-separated files you want to combine
## 2) filePattern --> the pattern in common to all files that will be used to grep them (eg, <your-first-sample>-summary.txt, <your-second-sample>-summary.txt )

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("ERROR: you need to specify the location of the directory containing the summary files for all your samples and their common ending filename ", call.=FALSE)
}

if(!require("dplyr")){
  r <- getOption("repos")
  r["CRAN"] <- "http://cran.us.r-project.org"
  options(repos = r)
  install.packages("dplyr", dependencies =TRUE)
}

inputDir = args[[1]]
filePattern = args[[2]]

if(endsWith(inputDir, '/')) {
  inputDir = paste(inputDir,'tables/',sep='')

} else {
  inputDir = paste(inputDir,'/tables/',sep='')
}

files = list.files(inputDir,recursive=T,full.names=T,pattern=filePattern)

fileNames = basename(files)

cat('\n Combining the below files \n')
fileNames
cat('located within: ', inputDir,'\n')


tables = lapply(files,function(f){
    f = read.table(f,header=T,sep='\t')
})

allTables <- do.call("rbind", tables)
outfile = paste(inputDir,'all-samples-summary-table.txt',sep='')
write.table(x = allTables,file = outfile,quote=F,row.names=F,col.names=T,sep='\t') 

cat('Combined table is located at: ',outfile ,'\n')



