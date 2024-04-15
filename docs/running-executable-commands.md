# Running executable commands

----------------------------------------------------------------

## Table of Contents
- [Running executable commands](#running-executable-commands)
  - [Table of Contents](#table-of-contents)
  - [Rationale](#rationale)
  - [Step-by-step guide](#step-by-step-guide)

----------------------------------------------------------------
## Rationale

Imagine of having some functions or tasks that need to be performed repeatedly. These can be, for instance, ensuring that the 2-weeks-no-touch file deletion policy of `/vast/scratch/` does not remove some files you need or retrieving all the archived files from a particular directory within `stornext` in order to (re)analyse some data. In all these situations, you have at least 3 choices:

1. Running the commands for each file and each time you need it (very time consuming and inconsistency/error-prone)
2. Creating a script in which you save the commands you want to execute (better) and only hard-code, ie, change, the input directory/file you want to process (improvable) 
3. Create a script and only pass command line arguments, ie flags, to it (easier)

As a general (bioinformatic) life principle: for any task that needs to be done more than once, just make have a script. Once you've create it, however, you might as well try to make your life even easier by soft-coding all the variables in a way that it will spare you from having to edit the script each time for every variable you need to change. This is not always possible, but for file management-related tasks (like the above examples) it can work pretty well. <br/>

For these reasons, I've created this document highlighting how to create an executable script and run it from any location within Milton. Having executables ready to use from anywhere, and soft-coded as much as possible so that you only need to specify command-line arguments, should greatly help you managing your own data and analyses. A collection of executables scripts can be found within the `bin/` folder of this repository. Below you can read a step-by-step guide on how to save and run any of those scripts on the HPC. <br/> 

----------------------------------------------------------------
## Step-by-step guide

1. Go into the `bin/` directory within this repository and identify the script you need (eg. `bin/resetTimestamp.sh`)
2. Click on that to open it and copy all the lines (eg, `CTRL +C`)
3. Log in into Milton from the command line (eg, `ssh vc7-shared` and enter your password)
4. You now will be in your `$HOME` directory (you can check its path by running `pwd`)
5. From here, if not already present, create a `bin` directory by executing  `mkdir $HOME/bin`
6. Open a text editor (eg, [nano](https://www.nano-editor.org/)) by running `nano $HOME/bin/resetTimestamp.sh`  
7. Paste in there all the lines you just copied (eg `CTRL +V`)
8. Save the results by typing `CTRL +O` and press `ENTER` then to exit the text editor type `CTRL +X`
9. Make the script executable by typing `chmod +x $HOME/bin/resetTimestamp.sh`. You'll need to run this command for each different script in the bin directory. For simplicity you can instead run `chmod +x $HOME/bin/*` which will make executable everything withing the `bin` directory. Remember this step otherwise you might get an error like this:`-bash: /home/users/allstaff/vespasiani.d/bin/resetTimestamp.sh: Permission denied`
10. Now that you have permessions to execute the script you need to make it executable from anywhere on the HPC. This means being able to call the script by simply typing `resetTimestamp.sh` instead of always having to specify its full location as `$HOME/bin/resetTimestamp.sh`. To this end, you need to export the path to your `$HOME/bin/` directory into your `bash_profile` which you can achieve by:
    * Opening up your `.bash_profile` configuration file by `nano $HOME/.bash_profile`
    * Exporting the path to the bin directory by writing this line in there `export PATH=$PATH:$HOME/bin/` &rarr this command will **permanently** add the new value of `PATH` to the shell environment so that each time you will login into the HPC and your bash profile is activated, your system will be able to locate that path
    * Saving and exiting nano text editor by typing `CTRL +O`, `ENTER` and then `CTRL +X`

To check that everything worked well you can run `command -v resetTimestamp.sh`. If everything went well, it will print the location of the file (which should be `$HOME/bin/`).
