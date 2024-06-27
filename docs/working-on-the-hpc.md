# Working on the HPC

----------------------------------------------------------------

## Table of Contents
- [Working on the HPC](#working-on-the-hpc)
  - [Table of Contents](#table-of-contents)
  - [Interactive sessions](#interactive-sessions)
  - [Modules](#modules)
  - [Batch jobs](#batch-jobs)
  - [Array jobs](#array-jobs)
  - [Resources](#resources)

----------------------------------------------------------------
## Interactive sessions

An interactive session is useful when you need to develop and test your code and/or perform quick data analyses. It allows you to interact with the HPC in real-time. <br/>
Interactive sessions are run on a computing node. So, from the log in node (ie, the computer you'll access after `ssh vc7-shared`), to start an interactive session you need to request computational resources. You can do it by running the following command from anywhere on the login node:

```
salloc --partition interactive --job-name "your-interactive-sess" --cpus-per-task 2  --time 24:00:00 --mem 50G
```

You can change any of those parameters depending on, eg:

* Which partition (aka subcluster of computers) you want to execute your job (check the "Slurm partition" post on the [Resources](#resources) section below)
* How much memory and time you think your job will need (sometimes this is a trial-error kinda of process) 
* Your job name
* etc..  

Note that each partition has its own limits in terms of time and memory available both in total and per user. So it can happen that your request won't be fulfilled because you are simply asking for too much or there are other users before you on the queue.

## Modules 

Modules are pre-installed softwares that you can load and use to run your jobs. Examples of modules can be R/python or suite of commands such as deeptools/samtools/bowtie2 etc... You can load a module by running  load a module you simply need to run `module load <modulename/module.version>` for example `module load deeptools/3.5.1`. Some modules may require others to be loaded first or as well. Most of the times this will be automatic, eg the R/4.4.0 module will automatically load the hdf5/1.8.21 module as well:
```
module load R/4.4.0
Loading R/4.4.0
  Loading requirement: hdf5/1.8.21
```
Other times this will not be the case and you'll figure it out when troubleshoot warnings and/or error messages from your script.  Because of this modules interdependency, module versions are extremely important to specify correctly and keep track. So while you can always load a module without specifying the software versions (eg `module load deeptools/`) and successfully run your analyses, you might not be able to do it anymore should newer software versions or dependencies will be installed on the HPC. Therefore always specify the version for each software you are using; otherwise be ready to change your script. <br/>

Note that installing softwares on the HPC is a prerogative of the system admin, ie WEHI IT support. If you need a new software to be installed for running your analyses you can ask IT to install it for you, although might be reluctant in doing it and likely suggest you other ways to install it. <br/>

To can check whether a module (and which version) is available on the HPC you can run `module avail`. This will return a list of all modules that can be loaded. Of course if you know the software you need is samtools, you can also run `module avail samtools` which will return all the possible versions of that module you can load,eg:
```
module avail samtools/ 
----------------------------------------------------------- /stornext/System/data/modulefiles/bioinf/its ------------------------------------------------------------
samtools/0.1.18  samtools/1.3.1  samtools/1.5  samtools/1.7  samtools/1.10  samtools/1.13  samtools/1.15    samtools/1.17  samtools/1.19.2  
samtools/0.1.19  samtools/1.4.1  samtools/1.6  samtools/1.9  samtools/1.12  samtools/1.14  samtools/1.16.1  samtools/1.18  samtools/1.20    
```

Once you have load the module (or a list of modules) you can check which one you have loaded by running `module list` and you can also unload a given module by running `module unload <modulename/module.version>` eg, `module unload deeptools/3.5.1`.

## Batch jobs
The other way to run jobs on the HPC is by submitting them to the scheduler (ie SLURM). SLURM is a managing system, ie assigns computational resources to users and manages job queues. Submitting jobs to the scheduler consists in specifying the above introduced `salloc` parameters into a bash script using a SLURM-specific syntax (see below).
```
!/bin/bash
#SBATCH --job-name=jobName 
#SBATCH --partition=regular
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.username@wehi.edu.au
#SBATCH --output=jobName-%A.out 
#SBATCH --error=jobName-%A.err

## Setup ------------------------------------------------------------------------
## Below you specify the list of modules you'd need to run the scripts and any eventual variable (eg, directory)
module load deeptools/3.5.1 
module samtools/1.20

## Your actual script ------------------------------------------------------------------------

## below you can specify each bash command you want to run.

```
You can [check this documentation](https://slurm.schedmd.com/sbatch.html) for a more complete list of all `#SBATCH` arguments that you can pass into the script (as well as to the `salloc` command).

Importantly, the above lines of code need to be entered and saved into a bash script (a file ending with `.sh`) which you can create from your terminal by using any text editor (eg, vim or nano). So after logging in into the HPC, you can open a text editor session by running `nano batch-script.sh` and copy the above lines of code into it. Then type `Ctrl + o` to save the edits into the `batch-script.sh` and `Ctrl + x` to exit the editing session.  

To execute your job you just need to run:
```
sbatch batch-script.sh
```
This command will submit your job to the SLURM scheduler which will allocate all the specified resourses, load all the requested modules and then run your script. Importantly, with this command you dont need to be constantly logged into the HPC. You can exit from the log in node and even shut your computer down and it will still be running for you. This because your script is now running (or is scheduled to be running) into a different computer (ie, a computing node).

Note that you can run `sbatch` from anywhere on the login node and also from a computing node (in case you are on an interactive session). What's important is that you correctly feed into `sbatch` the exact location (ie full path) to your `batch-script.sh`. Same applies for the lines of code within the actual script pointing to the input/output directories or files you want to use.

## Array jobs
(coming soon)

--------------------------------------------------
## Resources

Below you can find a list of resources you can consult to learn more about working on the HPC @ WEHI.

* [Slurm partitions](https://wehieduau.sharepoint.com/sites/rc2/SitePages/SLURM-partitions.aspx?web=1) 
* [Interactive workloads](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Interactive-workloads.aspx)
* [Batch jobs](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Getting-started-Slurm.aspx)
* [Array jobs](https://wehieduau.sharepoint.com/sites/rc2/SitePages/SLURM-job-arrays.aspx)
* [Advance job scripting](https://wehieduau.sharepoint.com/sites/rc2/SitePages/SLURM-advanced-job-scripts.aspx)