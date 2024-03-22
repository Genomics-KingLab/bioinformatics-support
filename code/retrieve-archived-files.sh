#!/bin/bash  


## this script will retrieve the disk copy of all archived files within an INPUT directory and its subdirectories
## copy the input directory into your $TMPDIR (ie, vast/scratch/users/$USER) 
## re-archive all files within directory and subdirectories in stornext

module load stornext/1.1

INPUTDIR="/stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/project_template" ## change this to whatever you want

function retrieveFilesIteratively() {
  local DIR="$1"

  for FILE in "$DIR"/*; do
    if [ -f "$FILE" ]; then

        echo
        echo "Retrieving archived $FILE"
        echo
        snretrieve "$FILE"

    fi

    if [ -d "$FILE" ]; then
      retrieveFilesIteratively "$FILE"
    fi
  done
}

retrieveFilesIteratively "$INPUTDIR"

echo
echo "No more files to retrieve"
echo "Copying everything into $TMPDIR"
echo

cp -r "$INPUTDIR" "$TMPDIR"


function rearchiveFilesIteratively() {
  local DIR="$1"

  for FILE in "$DIR"/*; do
    if [ -f "$FILE" ]; then
        echo
        echo "Re-archiving $FILE"
        echo
        snrmdiskcopy "$FILE"

    fi

    if [ -d "$FILE" ]; then
      rearchiveFilesIteratively "$FILE"
    fi
  done
}

rearchiveFilesIteratively "$INPUTDIR"
