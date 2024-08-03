#!/bin/bash  

#SBATCH --job-name=ENCODE-GM12878-WBS 
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=slurm-report/ENCODE-GM12878-WBS-%A.out
#SBATCH --error=slurm-report/ENCODE-GM12878-WBS-%A.err
#SBATCH --mem 200G 
#SBATCH --partition=bigmem

## Description: Script used to download ENCODE whole genome bisulfite sequencing processed data (bed file) for GM12878 cells and
## 1) Lift them to hs1 and 2) convert them to bigWigs for visualization

## Set up  ---------------------------------------------------------------------------------------------------------------------------------

module load ucsc-tools/
module load bedtools/

PROJECT="ENCODE_GM12878_WBS"
PROJECT_DIR="${TMPDIR}/$PROJECT"

if [ ! -d "$PROJECT_DIR" ]; then
  mkdir -p  "$PROJECT_DIR"
fi

cd $PROJECT_DIR

## get data and uncompress it
wget --content-disposition "https://www.encodeproject.org/files/ENCFF570TIL/@@download/ENCFF570TIL.bed.gz" && gunzip ENCFF570TIL.bed.gz 
wget --content-disposition "https://www.encodeproject.org/files/ENCFF614QHA/@@download/ENCFF614QHA.bed.gz" && gunzip ENCFF614QHA.bed.gz

## Convert bed to bedgraph  ---------------------------------------------------------------------------------------------------------------------------------

## convert bed to bedgraph with extra column being the % of methylation
cut -f 1,2,3,11 ENCFF570TIL.bed | grep -v "chrEBV" > GM12878_ENCFF570TIL_methylation.bed
cut -f 1,2,3,11 ENCFF614QHA.bed | grep -v "chrEBV" > GM12878_ENCFF614QHA_methylation.bed

# ## Liftover  ---------------------------------------------------------------------------------------------------------------------------------

## liftover the bedGraphs from hg38 to hs1
CHAIN="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/hg38/liftover/hg38ToHs1.over.chain"

SAMPLES=(
    "GM12878_ENCFF570TIL_methylation" 
    "GM12878_ENCFF614QHA_methylation"
)

## this takes a while
for s in "${SAMPLES[@]}";do
    liftOver -minMatch=0.95 -bedPlus=4 -tab ${s}.bed $CHAIN ${s}-hs1-mapped.bed ${s}-hs1-unmapped.bed
    sed '/^#/d' -i ${s}-hs1-unmapped.bed
done

## Convert bedgraph to bigwigs  ---------------------------------------------------------------------------------------------------------------------------------

CHROM_SIZES="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/hs1/annotations/hs1.chrom.sizes.sorted.txt"

## this also takes a whilex
for s in "${SAMPLES[@]}";do
    echo
    echo "Converting $s-hs1-mapped.bed from bedgraph to bigwig"
    echo
    sort -k1,1 -k2,2n $s-hs1-mapped.bed > $s-hs1-mapped-sorted.bed
    bedtools merge -i $s-hs1-mapped-sorted.bed -c 4 -d 0 -o max > $s-hs1-mapped-sorted-merged-max.bed ## merge overlapping regions and get their MAX value 
    bedtools merge -i $s-hs1-mapped-sorted.bed -c 4 -d 0 -o mean > $s-hs1-mapped-sorted-merged-mean.bed ## merge overlapping regions and get their MEAN value
    bedGraphToBigWig $s-hs1-mapped-sorted-merged-max.bed $CHROM_SIZES $s-hs1-mapped-sorted-merged-max.bigWig
    bedGraphToBigWig $s-hs1-mapped-sorted-merged-mean.bed $CHROM_SIZES $s-hs1-mapped-sorted-merged-mean.bigWig
done


echo
echo "Done"
echo