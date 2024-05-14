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
module load bowtie2


/stornext/General/data/academic/lab_king/BIOINFORMATICS/scripts/create-chimera-genome.sh -s ecoliASM584v2 -g hg38

