# Running executable commands

----------------------------------------------------------------

## Table of Contents
- [Running executable commands](#running-executable-commands)
  - [Table of Contents](#table-of-contents)
  - [Rationale](#rationale)
  - [Step by step guide](#step-by-step-guide)

----------------------------------------------------------------
## Rationale

Imagine of having some functions or tasks that need to be performed repeatedly. These can be, for instance, ensuring that the 2-weeks-no-touch file deletion policy of `/vast/scratch/` does not remove some **intermediate files** you are planning to use or retrieving all the archived files from a particular directory within `stornext` because you need to (re)analyse some data. In all these situations, you have at least 3 choices:

1. Running the commands for each file and each time you need it (very time consuming and inconsistency-prone)
2. Creating a script in which you save the commands you want to execute (better) and only hard-code, ie, change, the input directory/file you want to process (improvable) 
3. Create a script and only pass command line arguments, ie flags, to it (easier)

Because for any task that needs to be done more than once is best to have a script, once you create it you might as well try to make your life easier which means trying to avoid having to hard-code variables each time. This particularly works for file management related tasks, ie. all those actions you need to perform in order to best organise your files such as those examples mentioned above. For these reasons, I've created this document highlighting how to create an executable script and run it from any location within Milton. Also, the `bin/` folder of this repository contains a series of scripts that I've created that could be used by anyone within the King lab. Below you can read a step-by-step guide on how to save and run any of those scripts on the HPC. <br/> 

----------------------------------------------------------------
## Step by step guide

1. Go into the `bin/` directory within this repository and identify the script you need
2. Click on that to open it and copy all the lines (eg, `CTRL +C`)
3. Log in into Milton from the command line (eg, `ssh vc7-shared` and enter your password)
4. You now will be in your `$HOME` directory (you can check its path by running `pwd`)
5. From here create a bin directory by executing  `mkdir $HOME/bin`
6. Open a text editor (eg, [nano](https://www.nano-editor.org/)) by running `nano $HOME/bin/<the-same-filename-of-the-script-you-chose>`  
7. Paste in there all the lines you just copied (eg `CTRL +V`)
8. Save the results by typing `CTRL +O` and then exit the text editor by typing `CTRL +X`
9. Make the script executable by typing `chmod +x $HOME/bin/<the-same-filename-of-the-script-you-chose>`. You'll need to run this command for each different script in the bin directory. For simplicity you can instead run `chmod +x $HOME/bin/*` which will make executable everything withing the `bin` directory
10. Make the script executable from anywhere on Milton by exporting the path to the `$HOME/bin/` directory into your `bash_profile` by doing the following:
    * Open up your bash_profile configuration file by `nano $HOME/.bash_profile`
    * Export the path to the bin directory by writing this line in there `export PATH=$PATH:$HOME/bin/` &rarr this command will **permanently** add the new value of `PATH` to the shell environment so that each time you will login into the HPC and your bash profile is activated, your system will be able to locate that path
    * Saving and exiting nano text editor by typing `CTRL +O` and then `CTRL +X`

