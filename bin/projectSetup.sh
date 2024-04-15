
#!/bin/sh
## AUTHOR: Davide Vespasiani
## EMAIL: vespasiani.d@wehi.edu.au
## DESCRIPTION: use this script to set up all directories related to a new project

helpFunction()
{
   echo ""
   echo "Usage: "$0" -p projectName -sd stornextDir"
   echo -e "\t-p the name of your new project"
   echo -e "\t-sd path to the stornext dir where you want to create your project"
   exit 1 # Exit script after printing help
}

if (( $# == 0)); then
    echo "No input dir specified"
    echo "Remember to specify -i input dir with i being one of the following input directories";
    echo
    ls "$TMPDIR"
    echo
    echo
    echo "Exiting the script";
    exit 1
fi

while getopts "i:sd:" flag; do
    case "${flag}" in
        p) projectName=${OPTARG};;
        sd) stornextDir=${OPTARG};;
        ?) helpFunction ;; # Print helpFunction 
    esac
done


## make sure to be in the vast/scratch/$USER
CURRENTDIR=$(pwd)
TMPDIR="/vast/scratch/users/$USER"
INPUTDIR="$TMPDIR/$inputdir"


if [ ! -d "$INPUTDIR" ]; then

    echo
    echo "The specified $INPUTDIR directory does not exist."
    echo "Here is the list of all directories available in your tmp dir $TMPDIR"
    echo
    ls "$TMPDIR"
    echo "Exiting the script";
    exit 1; 

elif [[ "$CURRENTDIR" == "$TMPDIR" ]]; then
    echo
else 
    cd "$INPUTDIR"
fi

echo
echo "Resetting file timestamp for each file within the specified $INPUTDIR directory"
echo

find . -type f -exec touch {} +

echo
echo "Done, now your files have 2 weeks of life expectancy"
echo