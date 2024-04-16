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
    echo " -s = the path on stornext relative to /stornext/General/data/academic/lab_king/BIOINFORMATICS/$USER where you want to creat e your project"
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

if [[ $n != true && $sd != true ]]; then
    echo "Remember to specify both -n projectName -sd stornextDir flags in order to run this script"
    echo
    echo "Exiting the script";
    echo
    exit 1
fi 

## make sure to be in the vast/scratch/$USER
STORNEXT_PROJECTDIR="/stornext/General/data/academic/lab_king/BIOINFORMATICS/$USER/$projectName"
VASTSCRATCH_PROJECTDIR="/vast/scratch/users/$USER/$projectName"

OUTDIR='out'
DATADIR='data'
PLOTSDIR="$OUTDIR/plots"
FILESDIR="$OUTDIR/files"
TABLESDIR="$OUTDIR/tables"
CODEDIR='code'

directories=("$DATADIR" "$PLOTSDIR" "$FILESDIR" "$TABLESDIR" "$CODEDIR")

echo "Creating the project directories on stornext and on vast/scratch"
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