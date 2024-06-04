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
BASE_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/scripts/calibrated-pe-alignment" ## path to where calibrated-pe-alignment folder is located 
TMPDIR="/vast/scratch/users/$USER/" ## path to your vast/scrath tmp directory 
OUT_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/vespasiani.d/gabbyTest/calibrated-CutTag-Alignemnt-keepDup" ## path to where the output folder and files will be saved
DATA_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/vespasiani.d/gabbyTest/data/fastq" ## path to data directory containing the fastq files for your experiment
CHIMERAGENOME_DIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/databank/chimeras/ecoliASM584v2-hg38" ## path to the folder containing your chimera genome files
METADATAFILE="/stornext/General/data/academic/lab_king/BIOINFORMATICS/vespasiani.d/gabbyTest/data/metadata/gabbyTest-metadata.txt" ## path to the metadata txt file
METADATAFILENAME=$(basename ${METADATAFILE%.txt})
NEWMETADATAFILE="${TMPDIR}${METADATAFILENAME}_sample_${SLURM_ARRAY_TASK_ID}.txt"
SCRIPTS_DIR="${BASE_DIR}/scripts" 

echo
echo "Creating a temporary new metadata file by taking the line number $SLURM_ARRAY_TASK_ID of the original $METADATAFILENAME file and saving it at $NEWMETADATAFILE"
echo 

sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATAFILE > ${TMPDIR}${METADATAFILENAME}_sample_${SLURM_ARRAY_TASK_ID}.txt

echo
echo "Using the content of the new (single-line) metadata file to run the sbatch array job"
echo

${SCRIPTS_DIR}/calibrated-pe-alignment.sh -o" $OUT_DIR" -d "$DATA_DIR" -g "$CHIMERAGENOME_DIR"  -m "$NEWMETADATAFILE" -k yes

echo
echo "Combining all samples summary metadata information into a single file"
echo


Rscript ${SCRIPTS_DIR}/combine-summary-tables.R $OUT_DIR -summary.txt


echo
echo "Finished running the slurm job"
echo
