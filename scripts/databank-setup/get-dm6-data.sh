#!/bin/bash  

#SBATCH --job-name=getDm6data
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=getDm6data-%A.out
#SBATCH --error=getDm6data-%A.err
#SBATCH --mem 50G 

#Download genome files from UCSC
#Common name: human
#Taxonomic name: Homo sapiens, taxonomy ID: 9606
#Sequencing/Assembly provider ID: T2T Consortium
#Assembly date: 24 Jan 2022

module load pigz/
module load bowtie2/
module load samtools/
module load STAR/2.7.9a

GENOME='dm6'
DATABANK_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/$GENOME"

## Create directories  --------------------------------

if [ ! -d "$DATABANK_DIR/fasta" ]; then
    mkdir -p "$DATABANK_DIR/fasta"
fi

if [ ! -d "$DATABANK_DIR/genes" ]; then
    mkdir "$DATABANK_DIR/genes"
fi

if [ ! -d "$DATABANK_DIR/annotations" ]; then
    mkdir "$DATABANK_DIR/annotations"
fi

if [ ! -d "$DATABANK_DIR/2bit" ]; then
    mkdir "$DATABANK_DIR/2bit"
fi

if [ ! -d "$DATABANK_DIR/bowtie2" ]; then
    mkdir "$DATABANK_DIR/bowtie2"
fi

## Download files --------------------------------

## fasta
wget --content-disposition -P $DATABANK_DIR/fasta/ https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz

## 2bit index
wget --content-disposition -P $DATABANK_DIR/2bit/ https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.2bit 

## annotations
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chromAlias.bb
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chromAlias.txt

## genes
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ensGene.gtf.gz
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ncbiRefSeq.gtf.gz
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.refGene.gtf.gz

## ENCODE blacklisted regions
wget --content-disposition -P $DATABANK_DIR/annotations/ https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm6-blacklist.v2.bed.gz 

## Create extra files ---------------------------------------

## Chromosome sizes file
if (file $DATABANK_DIR/fasta/$GENOME.fa* | grep -q compressed ) ; then
    pigz --decompress $DATABANK_DIR/fasta/$GENOME.fa.gz
fi

## Bowtie2 index
bowtie2-build $DATABANK_DIR/fasta/$GENOME.fa $DATABANK_DIR/bowtie2/$GENOME

## Star index
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir $DATABANK_DIR/star \
--genomeFastaFiles $DATABANK_DIR/fasta/$GENOME.fa \
--sjdbGTFfile $DATABANK_DIR/genes/$GENOME.refGene.gtf
