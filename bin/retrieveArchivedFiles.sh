#!/bin/bash  
## AUTHOR: Davide Vespasiani
## EMAIL: vespasiani.d@wehi.edu.au
## DESCRIPTION: use this script to recursively retrieve the disk copy of all archived files within a given INPUT directory and its subdirectories
## copy the INPUT dir into your $TMPDIR (ie, vast/scratch/users/$USER) and then re-archive all files the INPUT directory 

module load stornext/1.1

KINGLAB_DIR="/stornext/General/data/academic/lab_king"
TMPDIR="/vast/scratch/users/$USER"

helpFunction()
{
    echo
    echo "Here's how you should run this script:  "$(basename $0)"  -i <inputDir> -o <outDir> -h "
    echo " -i = the path to the directory relative to $KINGLAB_DIR containing all the files you want to update"
    echo " -o = the path to the output directory relative to $TMPDIR where you want to save a copy of the files within the inputDir"
    echo " -h = to print this help function"
    echo   
    exit 1 # Exit script after printing help
}
a
if (( $# == 0)); then
    echo "No input dir specified"
    echo "Remember to specify -i input dir with i being the path relative to $KINGLAB_DIR to the subdirectory containing all the files you want to retrieve/archive";
    echo 
    echo "Exiting the script";
    exit 1
fi

while getopts ":hi:o:" flag; do
    case "${flag}" in
        h) helpFunction ;;
        i) inputdir=${OPTARG};;
        o) outdir=${OPTARG};;
       \?) echo "Invalid option: $OPTARG" >&2; exit 1 ;;
    esac
done

for FILE in $(find $inputdir -type f); do
    if [ -f "$FILE" ]; then
      echo
      echo "Retrieving $FILE"
      echo
      snretrieve "$FILE"
    fi
done

echo
echo "Copying the $inputdir directory and all its retrived files into the $outdir directory"
echo

cp -r $inputdir $outdir

echo
echo "Now re-archiving all files in the $inputdir directory"
echo

for FILE in $(find $inputdir -type f); do
    if [ -f "$FILE" ]; then
      echo
      echo "Archiving $FILE"
      echo
      snrmdiskcopy "$FILE"
    fi
done

echo
echo "All Done"
echo