#!/bin/bash  

#SBATCH --job-name=makeChimera
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=makeChimera.out  
#SBATCH --error=makeChimera.err  
#SBATCH --mem 50G
#SBATCH --time 1-00:00:00
#SBATCH --partition=regular

## check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi


## Modules -------------
module load bowtie2/2.4.4 
module load samtools/1.9

## specify I/O directories and files 
KINGLAB_DIR="/stornext/General/data/academic/lab_king"
SCRIPTS_DIR="${KINGLAB_DIR}/BIOINFORMATICS/pipes"  ## change this to where you have copied the align-to-chimera folder
SPIKEINGENOME_FASTA="${KINGLAB_DIR}/BIOINFORMATICS/databank/ecoliASM584v2/ecoliASM584v2.fa" ## point this to the fasta file for your spike-in genome
GENOMEOFINTEREST_FASTA="${KINGLAB_DIR}/BIOINFORMATICS/databank/hg38/hg38.fa"  ## point this to the fasta file for your genome of interest
OUT_DIR="${KINGLAB_DIR}/BIOINFORMATICS/databank/chimeras" ## change this to the location where you want to save the output of the script

${SCRIPTS_DIR}/align-to-chimera/scripts/create-chimera-genome.sh -s $SPIKEINGENOME_FASTA -g $GENOMEOFINTEREST_FASTA -o $OUT_DIR 

echo
echo "Finished running the job on slurm"
echo

