
#!/bin/sh
## AUTHOR: Davide Vespasiani
## EMAIL: vespasiani.d@wehi.edu.au
## DESCRIPTION: use this script to recursively touch each file in your vast/scratch directory in order to reset the timestamp

helpFunction()
{
   echo ""
   echo "Usage: $0 -i inputdir"
   echo -e "\t-i path to directory relative to vast/scratch/$USER containing files that need to be updated"
   exit 1 # Exit script after printing help
}

while getopts "i:" flag; do
    case "${flag}" in
        i) inputdir=${OPTARG};;
        ?) helpFunction ;; # Print helpFunction 
    esac
done

## make sure to be in the vast/scratch/$USER
CURRENTDIR=$(pwd)
TMPDIR="/vast/scratch/users/$USER"
INPUTDIR="$TMPDIR/$inputdir"

[[ -z "$i" ]] && { 
    echo
    echo "No input dir specified"
    echo "Remember to specify -i input dir";
    echo "Exiting the script";
    echo
    exit 1; 
}


if [[ "$CURRENTDIR" == "$TMPDIR" ]]; then
    echo
    echo
    echo
else 
    echo
    echo "Changing from current directory to $INPUTDIR"
    echo
    cd $INPUTDIR
fi

echo
echo "Resetting file timestamp for each file"
echo

find . -type f -exec touch {} +

echo
echo "Done, now your files have 2 weeks of life expectancy"
echo