#!/bin/bash  

#SBATCH --job-name=calibAlignkeepDup
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=slurm-report/calibAlignkeepDup-%A_%a.out
#SBATCH --error=slurm-report/calibAlignkeepDup-%A_%a.err
#SBATCH --mem 50G 

# Load modules: 
module load bowtie2/2.4.4 
module load samtools/1.9
module load R/4.4.0
module load bedtools/2.31.1
module load ucsc-tools/331

## specify I/O directories and files 
KINGLAB_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS"
TMPDIR="/vast/scratch/users/vespasiani.d/" ## path to your vast/scrath tmp directory
SCRIPTS_DIR="${KINGLAB_DIR}/scripts/calibrated-pe-alignment/scripts"  ## change this to where you have copied the calibrated-pe-alignment folder
OUT_DIR="${KINGLAB_DIR}/vespasiani.d/gabbyTest/calibrated-CutTag-Alignemnt-keepDup" ## change this to where you want the output files to be saved
DATA_DIR="${KINGLAB_DIR}/vespasiani.d/gabbyTest/data/fastq" ## change this to your data directory containing the fastq files
CHIMERAGENOME_DIR="${KINGLAB_DIR}/databank/chimeras/ecoliASM584v2-hg38" ## change this to the location of your chimera genome folder
METADATAFILE="${KINGLAB_DIR}/vespasiani.d/gabbyTest/data/metadata/gabbyTest-metadata.txt" ## This is a text file containing sample names 1 per row (and no header)

METADATAFILENAME=$(basename ${METADATAFILE%.txt})
NEWMETADATAFILE="${TMPDIR}${METADATAFILENAME}_sample_${SLURM_ARRAY_TASK_ID}.txt"

echo
echo "Creating a temporary new metadata file by taking the line number $SLURM_ARRAY_TASK_ID of the original $METADATAFILENAME file and saving it at $NEWMETADATAFILE"
echo 

sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATAFILE > ${TMPDIR}${METADATAFILENAME}_sample_${SLURM_ARRAY_TASK_ID}.txt

echo
echo "Using the content of the new (single-line) metadata file to run the sbatch array job"
echo

${SCRIPTS_DIR}/calibrated-pe-alignment.sh -o" $OUT_DIR" -d "$DATA_DIR" -g "$CHIMERAGENOME_DIR"  -m "$NEWMETADATAFILE" -k yes


echo
echo "Finished running the slurm job"
echo
