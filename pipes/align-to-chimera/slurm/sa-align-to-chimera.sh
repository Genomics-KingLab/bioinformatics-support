#!/bin/bash  

#SBATCH --job-name=align2chimera
#SBATCH --ntasks=1  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=slurm-report/align2chimera-%A_%a.out  
#SBATCH --error=slurm-report/align2chimera-%A_%a.err  
#SBATCH --mem 100G 

# Load modules: 
module load bowtie2/2.4.4 
module load samtools/1.9
module load R/4.4.0
module load bedtools/2.31.1
module load ucsc-tools/331

## specify I/O directories and files 
KINGLAB_DIR="/stornext/General/data/academic/lab_king"
SCRIPTS_DIR="${KINGLAB_DIR}/BIOINFORMATICS/pipes/align-to-chimera/scripts"
OUT_DIR="${KINGLAB_DIR}/BIOINFORMATICS/vespasiani.d/gabbyTest" ## change this to whathever you want
DATA_DIR="${KINGLAB_DIR}/RAW_DATA/gabby_tmp" ## point this to your data subdirectory within RAW_DATA/
CHIMERAGENOME_DIR="${KINGLAB_DIR}/BIOINFORMATICS/databank/chimeras/ecoliASM584v2-hg38" ## select your chimera genome of interest and remember to just add the file name (no '/' or file extensions)
METADATAFILE="${KINGLAB_DIR}/BIOINFORMATICS/vespasiani.d/gabbyTest/data/metadata/gabbyTest-metadata.txt" ## point this to your txt file containig the list of samples (1 per row with no header)
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $METADATAFILE) 

${SCRIPTS_DIR}/bowtie2-calibrated-pe-alignment.sh -o $OUT_DIR -d $DATA_DIR -s $SAMPLE -g $CHIMERAGENOME_DIR

echo
echo "Finished running the slurm array job"
echo

