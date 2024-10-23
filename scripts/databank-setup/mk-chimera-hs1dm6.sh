#!/bin/bash  

#SBATCH --job-name=hs1dm6-chimera
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=hs1dm6-chimera-%A.out
#SBATCH --error=hs1dm6-chimera-%A.err
#SBATCH --mem 50G 

module load pigz/
module load bowtie2/
module load samtools/
module load STAR/2.7.9a

ASSEMBLY_OI="hs1"
ASSEMBLY_SPIKE="dm6"
DATABANK_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes"
OUTDIR="${DATABANK_DIR}/${ASSEMBLY_OI}/chimeras"
LOGS_DIR="${OUTDIR}/logs"

DIR_GENOME_OI="${DATABANK_DIR}/${ASSEMBLY_OI}"
DIR_GENOME_SPIKE="${DATABANK_DIR}/${ASSEMBLY_SPIKE}"

## Outdir
CHIMERAGENOME="$ASSEMBLY_SPIKE-$ASSEMBLY_OI"
CHIMERAGENOME_DIR="${OUTDIR}/${CHIMERAGENOME}" 
LOGS_DIR="${OUTDIR}/logs" 

if [ ! -d "$LOGS_DIR" ]; then
  mkdir -p  "$LOGS_DIR"
fi

if [ ! -d "$CHIMERAGENOME_DIR" ]; then
  mkdir -p  "$CHIMERAGENOME_DIR"
fi

## Combine Fasta -----------------------------------------------------------------------------------------------------------------------------------

echo 
echo "Adding assembly identifier upstream of the chromosomes IDs contained within the $ASSEMBLY_SPIKE and the $ASSEMBLY_OI fasta files"
echo "And then concatenating the 2 fasta files into one"
echo

cp $DIR_GENOME_OI/fasta/$ASSEMBLY_OI.fa ${TMPDIR}/${ASSEMBLY_OI}.fa
cp $DIR_GENOME_SPIKE/fasta/$ASSEMBLY_SPIKE.fa ${TMPDIR}/${ASSEMBLY_SPIKE}.fa

## Add assembly name before each chromosome name in the fasta file
sed -i -e "s/>/>$ASSEMBLY_SPIKE-/g" "${TMPDIR}/${ASSEMBLY_SPIKE}.fa"
sed -i -e "s/>/>$ASSEMBLY_OI-/g" "${TMPDIR}/${ASSEMBLY_OI}.fa"

## create chimera fasta
cat ${TMPDIR}/${ASSEMBLY_OI}.fa ${TMPDIR}/${ASSEMBLY_SPIKE}.fa > ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa


## Index chimera fasta file -----------------------------------------------------------------------------------------------------------------------------------

echo 
echo "Indexing the chimera genome fasta file and creating a chrom.size.txt file from it"
echo

samtools faidx ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa -o ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa.fai
cut -f1,2 ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa.fai > ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.chrom.sizes.txt


## Combine GTF -----------------------------------------------------------------------------------------------------------------------------------
echo 
echo "Adding assembly identifier upstream of the chromosomes IDs contained within the $ASSEMBLY_SPIKE and the $ASSEMBLY_OI gtf files"
echo "And then concatenating the 2 gtf files into one"
echo

cp $DIR_GENOME_OI/genes/hs1.110.20220412.ncbiRefSeq.gtf ${TMPDIR}/hs1.110.20220412.ncbiRefSeq.gtf
cp $DIR_GENOME_SPIKE/genes/dm6.ncbiRefSeq.gtf.gz ${TMPDIR}/dm6.ncbiRefSeq.gtf.gz
gunzip ${TMPDIR}/dm6.ncbiRefSeq.gtf.gz

## Add assembly name before each gene name in the gtf file
sed -i "s/^/$ASSEMBLY_SPIKE-/" ${TMPDIR}/dm6.ncbiRefSeq.gtf
sed -i "s/^/$ASSEMBLY_OI-/" ${TMPDIR}/hs1.110.20220412.ncbiRefSeq.gtf

cat ${TMPDIR}/hs1.110.20220412.ncbiRefSeq.gtf ${TMPDIR}/dm6.ncbiRefSeq.gtf > ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.ncbiRefSeq.gtf

## Create Star index  -----------------------------------------------------------------------------------------------------------------------------------

STAR --runThreadN 30 \
--runMode genomeGenerate \
--genomeDir $CHIMERAGENOME_DIR/star \
--genomeFastaFiles $CHIMERAGENOME_DIR/$CHIMERAGENOME.fa \
--sjdbGTFfile $CHIMERAGENOME_DIR/$CHIMERAGENOME.ncbiRefSeq.gtf 

## Create Bowtie2 index  -----------------------------------------------------------------------------------------------------------------------------------

echo 
echo "Creating bowtie2 index for the $CHIMERAGENOME chimera genome"
echo

bowtie2-build  "${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa" "${CHIMERAGENOME_DIR}/${CHIMERAGENOME}" --threads 50 --verbose > "${LOGS_DIR}/bowtie2index-${CHIMERAGENOME}.log"

