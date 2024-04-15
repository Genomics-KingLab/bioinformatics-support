# Starting a computational project
----------------------------------------------------------------

## Table of Contents
- [Starting a computational project](#starting-a-computational-project)
  - [Table of Contents](#table-of-contents)
  - [Rationale](#rationale)

----------------------------------------------------------------

## Rationale
Best to spend one extra hour/day in organising the project before starting. Otherwise it will be very likely that you'll end up having soo many directories/subdirectories each containing figures, files,tables all over the places and with so many versions of each that you'll likely will not be able to find anything anymore. This will impact the quality, reproducibility and sharability of our research.
For this here are some guidelines to implement whenever starting a new research project.

1. Decide a project name, eg: `my-nobel-guaranteed-project`
2. Create a folder with this project name in each of the filesystems you have access to, eg. `stornext`, `vast/scratch` and, if it's a collaborative project, also on `vast/projects`. This means you should then have 2/3 project folders:

```
/vast/scratch/users/$USER/my-nobel-guaranteed-project  ## for all tmp files
/stornext/General/data/academic/lab_king/BIOINFORMATICS/$USER/my-nobel-guaranteed-project  ## for all final files and/or file versions that you want to lock
/vast/projects/my-nobel-guaranteed-project  ## for all files that will need to be shared/access by collaborators
```

Within each of these project folders you should have at least a `data/` and an `out/` subdirectory where to store (hard to believe) data and output files. 
Within the project folder located on stornext, you should also have at least a `code/` subdirectory where to store a copy of all your scripts (eg, your Rmarkdown files).
The `out/` subdirectory within each project folder should itself contain a `plots/`, `tables/` and `files/` directories in which, again, you'd save all the plots,tables and files, respectively. 
If you are starting to get confused and/ overwhelmed by the amount of housekeeping you need to do beforehand don't worry. Here is the overview of the ideal structure of each research project:

```
.
├── code
├── data
├── out
│   ├── files
│   ├── plots
│   └── tables

```
Also you can use the `bin/projectSetup.sh` script to create all these directories for each new project. Check out `docs/running-executable-commands.md` if you need to know how to create and run executables on the HPC.