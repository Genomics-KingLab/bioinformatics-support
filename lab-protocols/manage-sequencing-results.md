# Managing sequencing data
----------------------------------------------------------------
## Table of Contents
- [Managing sequencing data](#managing-sequencing-data)
  - [Table of Contents](#table-of-contents)
- [Where to start](#where-to-start)
  - [Workflow](#workflow)
  - [Examples](#examples)

----------------------------------------------------------------

# Where to start

So you sent some samples to be sequenced and got the results? Cool, then check the instructions and the workflow chart below as they contain information on how data is managed within the King lab. <br/>
If you need a refresher on the different file systems available on Milton and how they work you can check the documentation available at `docs/managing-files-and-disk-quota.md` still in this repository. <br>

Generally speaking, there are 2 main question you need to ask yourself before deciding where to store your files: What do you need to do with them and how desperate would you be if you'd loose them?
If those files are vital to you in that moment **AND** are hard/impossible to recover (eg you'd need to re-sequence samples), then save them on `stornext` and compress them. 
If you don't care about those files as you can recover them anytimes you want (eg, downloaded them from a public database), then just save them on `vast/scratch` (no need to compress them here). 
For any other intermediate situation you can save files on `vast/projects`, however this is not a unlimited storage nor in time nor in space.

## Workflow 
![Workflow file system management](../img/workflow-file-system.png)

## Examples

The table below gives you an overview on where to store the different files on the HPC depending on their type and/or their source.

| **Data**        	| **Type**     	| **Source**                                       	| **File example** 	| **Store**                     	| **Compress** 	| **Tar**    	|
|-----------------	|--------------	|--------------------------------------------------	|------------------	|-------------------------------	|--------------	|------------	|
| Public          	| Raw          	| Public repository accessible through script      	| Fastq            	| `vast/scratch`                	| No need      	| No need    	|
| Public          	| Raw          	| Public repository accessible through application 	| Fastq            	| `vast/project`                	| Yes          	| Preferably 	|
| Newly generated 	| Raw          	| Sequencing facilities                            	| Fastq            	| `storenext`                   	| Yes          	| Preferably 	|
| Intermediate    	| Preprocessed 	| Generate with custom scripts                     	| Bam              	| `vast/scratch`                	| No need      	| No need    	|
| Final           	| Preprocessed 	| Generate with custom scripts                     	| Bed, bigWig      	| `vast/project`<br>`storenext` 	| Yes          	| Preferably 	|
