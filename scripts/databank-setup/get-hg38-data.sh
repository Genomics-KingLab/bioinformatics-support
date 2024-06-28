#!/bin/bash  

#SBATCH --job-name=getHs1data
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=getHs1data-%A.out
#SBATCH --error=getHs1data-%A.err
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
module load cellranger/8.0.0

GENOME='hg38'
DATABANK_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/$GENOME"

## Create directories  --------------------------------

if [ ! -d "$DATABANK_DIR/fasta" ]; then
    mkdir "$DATABANK_DIR/fasta"
fi

if [ ! -d "$DATABANK_DIR/genes" ]; then
    mkdir "$DATABANK_DIR/genes"
fi

if [ ! -d "$DATABANK_DIR/liftover" ]; then
    mkdir "$DATABANK_DIR/liftover"
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
wget --content-disposition -P $DATABANK_DIR/fasta/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

## 2bit index
wget --content-disposition -P $DATABANK_DIR/2bit/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit 
wget --content-disposition -P $DATABANK_DIR/2bit/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit.bpt 

## annotations
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.bb
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes.txt

## genes
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

## liftover chains
wget --content-disposition -P $DATABANK_DIR/liftover/ https://hgdownload.soe.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz
wget --content-disposition -P $DATABANK_DIR/liftover/ https://hgdownload.soe.ucsc.edu/gbdb/hg38/liftOver/hg38ToHs1.over.chain.gz

## ENCODE blacklisted regions
wget --content-disposition -P $DATABANK_DIR/annotations/ https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz && ENCFF356LFX.bed.gz mv hg38_blacklist_ENCFF356LFX.bed.gz

## Create extra files ---------------------------------------

## Chromosome sizes file
if (file $DATABANK_DIR/fasta/$GENOME.fa* | grep -q compressed ) ; then
    pigz --decompress $DATABANK_DIR/fasta/$GENOME.fa.gz
fi

samtools faidx $DATABANK_DIR/fasta/$GENOME.fa | cut -f1,2 > $DATABANK_DIR/annotations/$GENOME.chrom.sizes.txt

## Bowtie2 index
bowtie2-build $DATABANK_DIR/fasta/$GENOME.fa $DATABANK_DIR/bowtie2/$GENOME

## Star index
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir $DATABANK_DIR/star \
--genomeFastaFiles $DATABANK_DIR/fasta/$GENOME.fa \
--sjdbGTFfile $DATABANK_DIR/genes/hg38.refGene.gtf


## Cellranger index
if (file $DATABANK_DIR/genes/hg38.refGene.gtf* | grep -q compressed ) ; then
    pigz --decompress $DATABANK_DIR/genes/hg38.refGene.gtf.gz
fi

cellranger mkref \
 --genome="$GENOME" \
 --fasta="${DATABANK_DIR}/${GENOME}.fa" \
 --genes="${DATABANK_DIR}/genes/hg38.refGene.gtf" \
 --nthreads=20 \
 --memgb=100 \
 --output-dir "${DATABANK_DIR}/cellranger"

