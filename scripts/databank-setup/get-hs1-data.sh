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

module load awscli/ 
module load pigz/
module load bowtie2/
module load samtools/
module load STAR/2.7.9a
module load cellranger/8.0.0

GENOME='hs1'
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
wget --content-disposition -P $DATABANK_DIR/fasta/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz 

## 2bit index
wget --content-disposition -P $DATABANK_DIR/2bit/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit 
wget --content-disposition -P $DATABANK_DIR/2bit/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit.bpt 

## annotations
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chromAlias.bb 
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chromAlias.txt 
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz 
wget --content-disposition -P $DATABANK_DIR/annotations/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.version.txt 

## genes
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gp.gz 
wget --content-disposition -P $DATABANK_DIR/genes/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz 

## files from https://github.com/marbl/CHM13?tab=readme-ov-file
aws s3 --no-sign-request cp s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz && mv chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz $DATABANK_DIR/genes/ ## JHU RefSeqv110 + Liftoff v5.1
aws s3 --no-sign-request cp s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3 && mv chm13.draft_v2.0.gene_annotation.gff3  $DATABANK_DIR/genes/ ## UCSC GENCODEv35 CAT/Liftoff v2

## liftover chains
wget --content-disposition -P $DATABANK_DIR/liftover/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg19.over.chain.gz
wget --content-disposition -P $DATABANK_DIR/liftover/ https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz 

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
--sjdbGTFfile $DATABANK_DIR/genes/chm13.draft_v2.0.gene_annotation.gff3 \
--sjdbGTFtagExonParentTranscript Parent ## this needs to be specified because i am using the gff3 file; generally default value is fine for gtf files 


## cellranger
## Note the gtf file used for this step is chm13.draft_v2.0.gene_annotation-gff3Converted-AttributesAndSourceFiltered.gtf
## The original file was the UCSC GENCODEv35 CAT/Liftoff v2 gff3 file (chm13.draft_v2.0.gene_annotation.gff3) which was:
## 1) Converted into a GTF using custom script (see bin/convert-gff2gtf.R)'
## 2) Filtered to remove multiple gene annotation sources (ie, CAT and Liftoff) --> in case of multiple sources I've kept only those with source = CAT; otherwise just the single reported annotation
## 2) Filtered to retain information for these attributes: 'protein_coding','lncRNA','IG_V_pseudogene','IG_V_gene','IG_D_gene','IG_J_gene','IG_J_pseudogene','IG_C_gene','IG_C_pseudogene','IG_pseudogene','TR_V_gene','TR_V_pseudogene','TR_J_gene','TR_J_pseudogene','TR_C_gene'

if (file $DATABANK_DIR/genes/chm13.draft_v2.0.gene_annotation-gff3Converted-AttributesAndSourceFiltered.gtf* | grep -q compressed ) ; then
    pigz --decompress $DATABANK_DIR/genes/chm13.draft_v2.0.gene_annotation-gff3Converted-AttributesAndSourceFiltered.gtf.gz
fi

cellranger mkref \
 --genome="$GENOME" \
 --fasta="${DATABANK_DIR}/${GENOME}.fa" \
 --genes="${DATABANK_DIR}/genes/chm13.draft_v2.0.gene_annotation-gff3Converted-AttributesAndSourceFiltered.gtf" \
 --nthreads=20 \
 --memgb=100 \
 --output-dir "${DATABANK_DIR}/cellranger"

