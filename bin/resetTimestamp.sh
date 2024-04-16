
#!/bin/sh
## AUTHOR: Davide Vespasiani
## EMAIL: vespasiani.d@wehi.edu.au
## DESCRIPTION: use this script to recursively touch each file in your vast/scratch directory in order to reset the timestamp and gain extra 2 weeks of life for them

helpFunction()
{
    echo
    echo "Here's how you should run this script:  "$(basename $0)"  -i inputdir -h "
    echo " -i = the path to the directory relative to /vast/scratch/$USER containing all the files you want to update"
    echo " -h = to print this help function"
    echo
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

while getopts ":hi:" flag; do
    case "${flag}" in
        h) helpFunction ;;
        i) inputdir=${OPTARG};;
       \?) echo "Invalid option: $OPTARG" >&2; exit 1 ;;
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