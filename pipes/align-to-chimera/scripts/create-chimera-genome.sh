#!/bin/sh

## Author: Hamish King and Davide Vespasiani
## email: king.h@wehi.edu.au and vespasiani.d@wehi.edu.au
## Description: Script to combine 2 genomes fasta files and create a: 1) fasta file (and its index); 2) bowtie2 index file;  3) chrom.size txt file for the resulting chimera and 4) chrom.size txt file for each of the 2 original genomes

helpFunction()
{
    echo  
    echo "Run this script as:  "$(basename $0)" -s <spike-in-genome> -g <genome-of-interest>  -o <outputDir> -h <helpFunction>"
    echo "-s = full path to the fasta file (including filename) of the spike-in genome. File can be gzip compressed and it must end with .fa.gz"
    echo "-g = full path to the fasta file (including filename) of the genome of interest. File can be gzip compressed and it must end with .fa.gz"
    echo "-o = the path to the directory where you want to save all the output of this script. Here the script will create 2 directories: 1) named <spike-in-genome>-<genome-of-interest> containing all output files and 2) a logs directory containing the log from bowtie2-build"
    echo "-h = to print this help function"
    echo
    exit 1 # Exit script after printing help
}

# Get the command line arguments
while getopts ":hs:g:o:" flag; do
    case "${flag}" in
        h) helpFunction ;;
        s) spikeGenome="${OPTARG}" ;; 
        g) genomeOfInterest="${OPTARG}" ;; 
        o) outputDir="${OPTARG}" ;; 
        \?) echo "Invalid option: ${OPTARG}" >&2; exit 1 ;;
    esac
done


## check if flags have been specified and if any of outputDir dir ends with / and, if so, remove the forward slash

if [[ -z "$spikeGenome" ]]; then
    echo
    echo "spikeGenome fasta file not specified"
    echo "Exiting the script";
    echo
    exit 1
fi 

if [[ -z "$genomeOfInterest" ]]; then
    echo
    echo "genomeOfInterest fasta file not specified"
    echo "Exiting the script";
    echo
    exit 1
fi 

if [[ -z "$outputDir" ]]; then
    echo
    echo "outputDir not specified"
    echo "Exiting the script";
    echo
    exit 1

 elif [[ "$outputDir" == */ ]]; then
      outputDir="${outputDir%/}"
fi 


## Set up dirs -----------------------------------------------------------------------------------------------------------------------------------
SCRIPT_ID="creating-chimera-genome"
SPIKEINGENOME_FASTA="$spikeGenome"
GENOMEOFINTEREST_FASTA="$genomeOfInterest"
TMP_DIR="/vast/scratch/users/${USER}/$SCRIPT_ID"
LOGS_DIR="${outputDir}/logs/"


if [ ! -d "$TMP_DIR" ]; then
  mkdir -p  "$TMP_DIR"
fi

if [ ! -d "$LOGS_DIR" ]; then
  mkdir -p  "$LOGS_DIR"
fi


if [[ $SPIKEINGENOME_FASTA == *.fa.gz ]]; then

    SPIKEINGENOME="$(basename "${SPIKEINGENOME_FASTA%%.fa.gz}")"
    gunzip -c "$SPIKEINGENOME_FASTA" > "${TMP_DIR}/${SPIKEINGENOME}.fa"
  else
    SPIKEINGENOME="$(basename "${SPIKEINGENOME_FASTA%.fa}")"
    cp "$SPIKEINGENOME_FASTA" "${TMP_DIR}/${SPIKEINGENOME}.fa"

fi

if [[ $GENOMEOFINTEREST_FASTA == *.fa.gz ]]; then

    GENOMEOFINTEREST="$(basename "${GENOMEOFINTEREST_FASTA%%.fa.gz}")"
    gunzip -c "$GENOMEOFINTEREST_FASTA" > "${TMP_DIR}/${GENOMEOFINTEREST}.fa"
  else
    GENOMEOFINTEREST="$(basename "${GENOMEOFINTEREST_FASTA%.fa}")"
    cp "$GENOMEOFINTEREST_FASTA" "${TMP_DIR}/${GENOMEOFINTEREST}.fa"

fi


CHIMERAGENOME="$SPIKEINGENOME-$GENOMEOFINTEREST"
CHIMERAGENOME_DIR="${outputDir}/${CHIMERAGENOME}"

if [ ! -d "$CHIMERAGENOME_DIR" ]; then
  mkdir -p  "$CHIMERAGENOME_DIR"
fi

##------------------------------------------------------------------------------------------------------------------------------------

echo 
echo "Adding assembly identifier upstream of the chromosomes IDs contained within the $SPIKEINGENOME and the $GENOMEOFINTEREST fasta files"
echo

sed -i -e "s/>/>$SPIKEINGENOME-/g" "${TMP_DIR}/${SPIKEINGENOME}.fa"
sed -i -e "s/>/>$GENOMEOFINTEREST-/g" "${TMP_DIR}/${GENOMEOFINTEREST}.fa"

##------------------------------------------------------------------------------------------------------------------------------------

echo 
echo "Concatenating the $SPIKEINGENOME and the $GENOMEOFINTEREST genomes"
echo

cat ${TMP_DIR}/${GENOMEOFINTEREST}.fa ${TMP_DIR}/${SPIKEINGENOME}.fa > ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa

##------------------------------------------------------------------------------------------------------------------------------------

echo 
echo "Creating bowtie2 index for the $CHIMERAGENOME chimera genome"
echo

bowtie2-build  "${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa" "${CHIMERAGENOME_DIR}/${CHIMERAGENOME}" --threads 50 --verbose > "${LOGS_DIR}/bowtie2index-${CHIMERAGENOME}.log"

##------------------------------------------------------------------------------------------------------------------------------------

echo 
echo "Indexing the chimera genome fasta file and creating a chrom.size.txt file from it"
echo

samtools faidx ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa -o ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa.fai
cut -f1,2 ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.fa.fai > ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.chrom.sizes.txt

##------------------------------------------------------------------------------------------------------------------------------------

echo 
echo "Extracting the chromosome sizes also for each of the 2 genomes separately"
echo

grep $GENOMEOFINTEREST ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.chrom.sizes.txt | sed -r s/$GENOMEOFINTEREST-//g  > ${CHIMERAGENOME_DIR}/${GENOMEOFINTEREST}.chrom.sizes.txt
grep $SPIKEINGENOME ${CHIMERAGENOME_DIR}/${CHIMERAGENOME}.chrom.sizes.txt | sed -r s/$SPIKEINGENOME-//g  > ${CHIMERAGENOME_DIR}/${SPIKEINGENOME}.chrom.sizes.txt

##------------------------------------------------------------------------------------------------------------------------------------

echo 
echo "All done!"
echo
echo "All output files are in the $CHIMERAGENOME_DIR"
echo
