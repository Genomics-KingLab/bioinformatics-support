#!/bin/bash  

#SBATCH --job-name=align2chimera
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=align2chimera-%A_%a.out  
#SBATCH --error=align2chimera-%A_%a.err  
#SBATCH --mem 50G 

# Load modules: 
module load bowtie2/2.4.4 
module load samtools/1.9
module load R/4.4.0
module load bedtools/2.31.1
module load ucsc-tools/331

## specify I/O directories and files 
KINGLAB_DIR="/stornext/General/data/academic/lab_king"
SCRIPTS_DIR="${KINGLAB_DIR}/BIOINFORMATICS/pipes"  ## change this to where you have copied the align-to-chimera folder
OUT_DIR="${KINGLAB_DIR}/BIOINFORMATICS/vespasiani.d/gabbyTest" ## change this to where you want the output files to be saved
DATA_DIR="${KINGLAB_DIR}/RAW_DATA/gabby_tmp" ## change this to your data directory containing the fastq files
CHIMERAGENOME_DIR="${KINGLAB_DIR}/BIOINFORMATICS/databank/chimeras/ecoliASM584v2-hg38" ## change this to the location of your chimera genome folder
METADATAFILE="${KINGLAB_DIR}/BIOINFORMATICS/vespasiani.d/gabbyTest/data/metadata/gabbyTest-metadata.txt" ## This is a text file containing sample names 1 per row (and no header)
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $METADATAFILE) 

${SCRIPTS_DIR}/align-to-chimera/scripts/bowtie2-calibrated-pe-alignment.sh -o $OUT_DIR -d $DATA_DIR -s $SAMPLE -g $CHIMERAGENOME_DIR

echo
echo "Finished running the slurm job"
echo

