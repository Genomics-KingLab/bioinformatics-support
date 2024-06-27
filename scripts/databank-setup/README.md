# Content
This folder contains bash scripts wrapping the list of commands used to download the data for the king-lab data bank directory on stornext.
Each script is meant to be executed as a sbatch job on slurm given steps such as bowtie2/star/cellranger indexes can take up long time. These scripts are also meant to be executed only once unless major changes happen at source. 

