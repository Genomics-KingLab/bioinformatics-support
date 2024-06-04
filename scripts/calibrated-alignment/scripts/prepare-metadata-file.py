#!/usr/bin/env python

## Author: Hamish King and Davide Vespasiani
## email: king.h@wehi.edu.au and vespasiani.d@wehi.edu.au
## Description: This script receives a csv (or xlsx) file prepares a tab-separated metadata file in a format acceptable for calibrated-pe-alignment.sh to run
## Usage: Call this script from the command line specifying as arguments (in order):
## 1) input --> the full path to the input file 
## 2) Sample column name
## 3) Read 1 column name
## 4) Read 2 column name

import os 
import pandas as pd
import sys

inputFile = sys.argv[1]
sampleCol = sys.argv[2]
read1Col = sys.argv[3]
read2Col = sys.argv[4]

if inputFile.endswith('.xlsx'):
   metadata = pd.read_excel(inputFile, sheet_name=0)
else:
   metadata = pd.read_csv(inputFile)

metadata = metadata[[sampleCol, read1Col, read2Col]]

outdir = os.path.dirname(inputFile)
outfileName = outdir + '/' + os.path.splitext(os.path.basename(inputFile))[0] + ".txt"

print("\n Exporting the metadata file to {} \n".format(outfileName))
print("\n Specify this when running the calibrated-pe-alignment.sh script \n")
metadata.to_csv(outfileName,sep='\t',header=None,index=None)





