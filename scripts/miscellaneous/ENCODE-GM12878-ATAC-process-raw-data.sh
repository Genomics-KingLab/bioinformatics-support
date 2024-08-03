#!/bin/bash  

#SBATCH --job-name=ENCODE-GM12878-ATAC 
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=slurm-report/ENCODE-GM12878-ATAC-%A_%a.out
#SBATCH --error=slurm-report/ENCODE-GM12878-ATAC-%A_%a.err
#SBATCH --mem 200G 
#SBATCH --partition=bigmem
#SBATCH --array=1-2

## Description: Script to download and process the ATAC-seq data related to the ENCODE experiment ENCSR095QNB (https://www.encodeproject.org/experiments/ENCSR095QNB/).
## Experiment details: ATAC-seq data for 2 biological replicates (ENCLB907YRF and ENCLB584REF) of GM12878 cell line


## Set up  ---------------------------------------------------------------------------------------------------------------------------------

module load sra-toolkit/ 
module load pigz/
module load samtools/
module load bowtie2/ 
module load deeptools/3.5.5
module load bedtools/

PROJECT="ENCODE_GM12878_ATAC"
INPUTFILE="/vast/scratch/users/vespasiani.d/$PROJECT/array-files/sra-list.txt" ## change this to your input file (2 columns tab-separated file no header). col1 = SRA ID, col2 = sample name
TMPDIR="/vast/scratch/users/vespasiani.d" ## change this to your own directory on vast/scratch
PROJECT_DIR="${TMPDIR}/$PROJECT"
DATA_DIR="${PROJECT_DIR}/data"
OUT_DIR="${PROJECT_DIR}/out"
LOG_DIR="${PROJECT_DIR}/log"
BAM_DIR="${OUT_DIR}/bam"
BED_DIR="${OUT_DIR}/bed"
BIGWIG_DIR="${OUT_DIR}/bigwig"

DIRS=("$PROJECT_DIR" "$DATA_DIR" "$OUT_DIR" "$LOG_DIR" "$BAM_DIR" "$BED_DIR" "$BIGWIG_DIR")

for dir in "${DIRS[@]}"; do
  if [ ! -d "$dir" ]; then
    mkdir  "$dir"
  fi
done

cd $PROJECT_DIR

SRA="$(sed -n "$SLURM_ARRAY_TASK_ID"p $INPUTFILE | awk '{print $1}' )" ## takes first column of line(s) = 1-n of the job array
SAMPLE="$(sed -n "$SLURM_ARRAY_TASK_ID"p $INPUTFILE | awk '{print $2}' )" ## takes second column of line(s) = 1-n of the job array

## Download raw data ---------------------------------------------------------------------------------------------------------------------------------

R1="${DATA_DIR}/${SRA}_1.fastq"
R2="${DATA_DIR}/${SRA}_2.fastq"

if [[ ! -f $R1 && -f $R2 ]]; then
  echo
  echo "Downloading fastq files for $SRA"
  echo
  fastq-dump.3.0.0 --outdir $DATA_DIR  --split-files $SRA
fi

## Align fastq reads -----------------------------------------------------------------------------------------------------------------------------------

ASSEMBLY="hs1"
BOWTIE2_INDEX="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/$ASSEMBLY/bowtie2/$ASSEMBLY"

echo
echo "Aligning samples to $ASSEMBLY genome assembly"
echo

bowtie2 -p 100 -q --sensitive --no-mixed --no-discordant --dovetail --mm -x $BOWTIE2_INDEX -1 $R1 -2 $R2  2> ${LOG_DIR}/${SAMPLE}-bowtie2-alignment-results.txt | grep -v XS: - | samtools view -bhS -F4 - > ${BAM_DIR}/${SAMPLE}-reads-uniqMap.bam

## Filter aligned reads  -----------------------------------------------------------------------------------------------------------------------------------

echo
echo "Removing duplicated reads and those mapping to chrM"
echo

samtools view -h ${BAM_DIR}/${SAMPLE}-reads-uniqMap.bam | grep -v 'chrM' | samtools sort -n -o ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM.bam

## sort bams by name (required for fixmate)
samtools sort --threads 20 -n ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM.bam -o ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName.bam  
samtools fixmate --threads 20 -m ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName.bam  ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName-fixmate.bam 

## resort bams by coordinates (required for identifying and removing duplicates)
samtools sort --threads 20 ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName-fixmate.bam  -o ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortPos-fixmate.bam 

## remove duplicates
samtools markdup --threads 20 -r ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortPos-fixmate.bam  ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortPos-fixmate-dedup.bam

## resort bams by name
samtools sort --threads 20 -n ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortPos-fixmate-dedup.bam -o ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName-fixmate-dedup-final.bam ## file to keep

## Convert PE reads to fragments -----------------------------------------------------------------------------------------------------------------------------------

echo
echo "Creating a fragment bed file from the uniquely mapped and deduplicated reads"
echo

bedtools bamtobed -i ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName-fixmate-dedup-final.bam -bedpe > ${BED_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName-fixmate-dedup-final.bed

## This command takes the 1st, 2nd and 6th column of the bedpe file (respectively containing the chr name, 5' min start and 3' max end of the 4 reads) and creates a fragment file with the resulting coordinates
cut -f 1,2,6,7 ${BED_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName-fixmate-dedup-final.bed | sort -k1,1 -k2,2n -k3,3n > ${BED_DIR}/${SAMPLE}-fragments-uniqMap-noChrM-sortName-fixmate-dedup-final.bed 

echo
echo "Converting the fragments bed file to a fragments bam file"
echo

CHROM_SIZES="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/genomes/hs1/annotations/hs1.chrom.sizes.sorted.txt"

bedToBam -i ${BED_DIR}/${SAMPLE}-fragments-uniqMap-noChrM-sortName-fixmate-dedup-final.bed  -g $CHROM_SIZES > ${BAM_DIR}/${SAMPLE}-fragments-uniqMap-noChrM-sortName-fixmate-dedup-final.bam ## file to keep

## Create genome coverage tracks -----------------------------------------------------------------------------------------------------------------------------------

echo
echo "Creating a genome coverage bedGraph file for the fragments that uniquely mapped to hs1"
echo

## Calculating the genome coverage without scaling the number of reads
genomeCoverageBed -bga -ibam ${BAM_DIR}/${SAMPLE}-fragments-uniqMap-noChrM-sortName-fixmate-dedup-final.bam | sort -k1,1 -k2,2n  > ${BED_DIR}/${SAMPLE}-fragments-uniqMap-noChrM-sortName-fixmate-dedup-final-genomeCoverage.bedGraph

echo
echo "Converting bedGraph to bigWig"
echo

bedGraphToBigWig ${BED_DIR}/${SAMPLE}-fragments-uniqMap-noChrM-sortName-fixmate-dedup-final-genomeCoverage.bedGraph $CHROM_SIZES ${BIGWIG_DIR}/${SAMPLE}-fragments-uniqMap-noChrM-sortName-fixmate-dedup-final-genomeCoverage.bigWig ## file to keep


## Remove unnecessary files -----------------------------------------------------------------------------------------------------------------------------------

rm ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM.bam
rm ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName.bam 
rm ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortName-fixmate.bam 
rm ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortPos-fixmate.bam 
rm ${BAM_DIR}/${SAMPLE}-reads-uniqMap-noChrM-sortPos-fixmate-dedup.bam
rm ${BED_DIR}/* ## no need to keep any of the bed files as you already have bams and bigwig

