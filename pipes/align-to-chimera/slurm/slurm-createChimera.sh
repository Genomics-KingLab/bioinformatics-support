#!/bin/bash  

#SBATCH --job-name=createChimera
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=slurm-report/chimera-genome-%A_%a.out  
#SBATCH --error=slurm-report/chimera-genome-%A_%a.err  
#SBATCH --mem 100G
#SBATCH --time 2-00:00:00
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
SPIKEINGENOME_FASTA="${KINGLAB_DIR}/BIOINFORMATICS/databank/ecoliASM584v2/ecoliASM584v2.fa" ## change this to your data directory containing the fastq files
GENOMEOFINTEREST_FASTA="${KINGLAB_DIR}/BIOINFORMATICS/databank/hg38/hg38.fa" ## change this to your data directory containing the fastq files
OUT_DIR="${KINGLAB_DIR}/BIOINFORMATICS/databank/chimeras" ## change this to the location of your chimera genome folder

${SCRIPTS_DIR}/align-to-chimera/scripts/create-chimera-genome.sh -s $SPIKEINGENOME_FASTA -g $GENOMEOFINTEREST_FASTA -o $OUT_DIR 

echo
echo "Finished running the job on slurm"
echo

