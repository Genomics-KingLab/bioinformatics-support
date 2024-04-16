#!/bin/sh
## AUTHOR: Davide Vespasiani
## EMAIL: vespasiani.d@wehi.edu.au
## DESCRIPTION: use this script to set up all directories related to a new project

# Define the help function
helpFunction()
{
    echo
    echo "Here's how you should run this script:  "$(basename $0)" -n projectName -s stornextDir -h "
    echo " -n = the name of your new project"
    echo " -s = the path to stornext where you want to create the project, by default this is /stornext/General/data/academic/lab_king/BIOINFORMATICS/$USER/"
    echo " -h = to print this help function"
    echo
    exit 1 # Exit script after printing help
}


# Get the command line arguments
while getopts ":hn:s:" flag; do
    case "${flag}" in
        h) helpFunction ;; 
        n) projectName="${OPTARG}" ;;
        s) stornextDir="${OPTARG}" ;;
        \?) echo "Invalid option: $OPTARG" >&2; exit 1 ;;
    esac
done

if [[ -z "$projectName" ]]; then
    echo
    echo "Project name not specified"
    echo "To run this script remember to define at least the projectName by specifying the -n flag"
    echo "Exiting the script";
    echo
    exit 1
fi 

## make sure to be in the vast/scratch/$USER

if [[ ! -z "$stornextDir" ]]; then

    STORNEXT_PROJECTDIR="$stornextDir"
else 
    STORNEXT_PROJECTDIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/$USER/$projectName"

fi 

VASTSCRATCH_PROJECTDIR="/vast/scratch/users/$USER/$projectName"

OUTDIR='out'
DATADIR='data'
PLOTSDIR="$OUTDIR/plots"
FILESDIR="$OUTDIR/files"
TABLESDIR="$OUTDIR/tables"
CODEDIR='code'

directories=("$DATADIR" "$PLOTSDIR" "$FILESDIR" "$TABLESDIR" "$CODEDIR")

echo "Creating the project directories on: "
echo " 1. stornext -> $STORNEXT_PROJECTDIR"
echo " 2. vast/scratch -> $VASTSCRATCH_PROJECTDIR"
echo

for d in "${directories[@]}"; do

  mkdir -p "$STORNEXT_PROJECTDIR/$d"
  mkdir -p "$VASTSCRATCH_PROJECTDIR/$d"

done

if [[ $c == true ]]; then
    
    for d in "${directories[@]}"; do
        mkdir -p "$STORNEXT_PROJECTDIR/$d"
        mkdir -p "$VASTSCRATCH_PROJECTDIR/$d"
    done
fi 


echo
echo "Done, your project folders are set up"
echo