# Starting a computational project
----------------------------------------------------------------

## Table of Contents
- [Starting a computational project](#starting-a-computational-project)
  - [Table of Contents](#table-of-contents)
  - [Rationale](#rationale)
  - [Set up a new project on Milton](#set-up-a-new-project-on-milton)

----------------------------------------------------------------

## Rationale

When it comes to data analysis, from experience, it is best to spend extra time in structuring your project folder before starting. Otherwise it is highly likely that you'll end up having soo many directories/subdirectories or (worse) figures, files, tables all over the places. Because bioinformatics (and science) is an iterative profession, not having an organised structure will result in significant time wasted on retrieving the right file, with the strong chance of inconsistency and having to repeat everything from scratch. Organising a project folder before starting will therefore improve the quality, reproducibility and sharability of our research but, most importantly, it will save you soo much time and trouble.

For all these reasons, here you can find some guidelines to implement whenever starting a new research project within the lab (and beyond). <br/>
**Note**, these do not represent the only way things have to be done. If you want to organize your research differently, you should do it. However, quoting Nike: "JUST DO IT".

1. Decide a project name, eg: `my-nobel-guaranteed-project`
2. Create a folder with this project name in each of the filesystems you have access to, eg. `stornext`, `vast/scratch` and, if it's a collaborative project, also on `vast/projects`. This means you should then have 2/3 project folders:

```
/vast/scratch/users/$USER/my-nobel-guaranteed-project  ## for all tmp files
/stornext/General/data/academic/lab_king/BIOINFORMATICS/$USER/my-nobel-guaranteed-project  ## for all final files and/or file versions that you want to lock
/vast/projects/my-nobel-guaranteed-project  ## for all files that will need to be shared/access by collaborators
```

Within each of these project folders you should have at least a `data/` , a `code/` and an `out/` sub-directories. As you can imagine:
1. `data/` &arr the folder where you would store all the input files
2. `code/` &arr the folder where you would store all your scripts (eg, your Rmarkdown files)
3. `out/` &arr the folder where you would store all the outputs of your analyses
The `out/` subdirectory within each project folder should itself contain a `plots/`, `tables/` and `files/` sub-directories where, again, you'd save all the plots, tables and files, respectively. 

Here is the overview of an ideal structure each research project folder should have:

```
.
├── code
├── data
├── out
│   ├── files
│   ├── plots
│   └── tables

```
----------------------------------------------------------------

## Set up a new project on Milton

To set up a new project-folder with the above pictured structure you can use the `bin/projectSetup.sh` script. Copy it into your `$HOME/bin` directory on Milton (if you don't know how, read the `docs/running-executable-commands.md` guidelines). Once copied and ensured it can run, you can create as many project folders as you like by running `projectSetup.sh -n my-nobel-guaranteed-project`.  The script will create a project directory on vast/scratch (ie `/vast/scratch/$USER/my-nobel-guaranteed-project`) and on stornext (ie, `/stornext/General/data/academic/lab_king/BIOINFORMATICS/$USER/my-nobel-guaranteed-project`), both containing the above subdirectories. 
