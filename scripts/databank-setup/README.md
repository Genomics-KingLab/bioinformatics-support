# Content

This folder contains all scripts used to populate the databank directory within the the King-lab space on stornext.
These scripts either get the data from UCSC Genome Browser (or other online databases) or use the downloaded data to create chimeric genomes which we need for calibrated alignments Each script is meant to be executed as a sbatch job on slurm given steps such as bowtie2/star/cellranger indexes can take up long time. 
Also, these scripts are meant to be executed only once unless major changes happen at source.

