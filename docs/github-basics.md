# Github basics
----------------------------------------------------------------
## Table of Contents
- [Github basics](#github-basics)
  - [Table of Contents](#table-of-contents)
  - [Overview and scope](#overview-and-scope)
  - [Configuration](#configuration)
    - [Setting up SSH key](#setting-up-ssh-key)
  - [Basic git usage](#basic-git-usage)

----------------------------------------------------------------

## Overview and scope

For a detailed information on github and its basic usage you can check [this post](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Basic-Git-Usage.aspx) and [this repo](https://github.com/WEHI-ResearchComputing/GitIntro/blob/master/Lessons/GitIntro_1.md) from WEHI IT support. In a nutshell, github allows you to create, manage and share computation notebooks. This benefits you manifolds as you can:

1. Safely store and create a back up for your code 
2. Keep track of each change 
3. Share your code with other people, thus improving collaborations, reproducibility and openess of science

## Configuration

Assuming you have successfully logged into github (otherwise you shouldnt be able to read this document), below you can read a step-by-step guide on how to configure your laptop to allow you to create repositories, push and pull changes into the KingLab github account.

### Setting up SSH key

All the information reported here on how to set up github with ssh keys come from [this post](https://gist.github.com/xirixiz/b6b0c6f4917ce17a90e00f9b60566278) and [this github page](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent). <br/>

1. Create a new ssh key to be linked with the github email address by running this command either on Windows or macOS: `ssh-keygen -t rsa -b 4096 -C "genomicskinglab@gmail.com"` and press `Enter` when the terminal will prompt the following:

```
Enter file in which to save the key (~/.ssh/id_rsa):
Enter passphrase (empty for no passphrase):  
Enter same passphrase again:
```

2. Copy the newly generated key on your keyboad by running `pbcopy < ~/.ssh/id_rsa.pub` on macOS or `cat ~/.ssh/id_rsa.pub | clip` on Windows. **NB:** this is just like running `CTRL + C` on a selected text. Then go into the lab's github.com account > Settings > SSH and GPG keys, click onto the `New SSH key` button, copy the key in there (ie `CTRL + V`) and give it a meaningful title (eg `<your-name>-ssh-authentication-key`).

3. Go back to your terminal and run `ssh -T git@github.com` to check whether the ssh key successfully works and you can be authenticated into the lab's github account. The terminal should prompt you something like this:

```
Hi Genomics-KingLab! You've successfully authenticated, but GitHub does not provide shell access.
```

4. Add the ssh key to the agent by running `ssh-add ~/.ssh/id_rsa`. **NB:** this might not work on Windows (not sure why).

5. Once you successfully connected your local machine to github you can start adding your folders in the Kinglab github account. See the next section for more information on how to use github.
   
## Basic git usage

The steps below represent the very basic commands you need to use in order to syncronise any local folder with the github account. <br/>

Imagine having/creating a folder on your local computer named `myProject` with the below structure:
```
$ myProject
.
├── code
│   ├── script1.R
│   └── script2.py
├── data
│   ├── file1.fastq
│   └── file2.bam
└── README.md
```

In order to syncronise this project to the github account and start keeping track of all the changes happening to the files within the directory, you need to:

* `cd myProject` into your project folder from the terminal
* Initialise the `myProject` directory as a git repository by running `git init`
* Go into the [KingLab github account](https://github.com/Genomics-KingLab) then > Repositories > New
* Give the new repository a name (eg again `myProject`) and make it private (you can change it to public anytime)
* Follow the command line instructions that appear. If you arleady have initialised (`git init`) a repository on your local computer and you just want to push it to github, then type the following commands on the terminal (make sure to be in your `myProject` directory!):

```
git remote add origin git@github.com:Genomics-KingLab/myProject.git
git branch -M main
git push -u origin main
```

* Check that your local folder is linked to the online github repository by running `git remote -v` which should print:

```
origin  git@github.com:Genomics-KingLab/myProject.git (fetch)
origin  git@github.com:Genomics-KingLab/myProject.git (push)
```

* Once this is done, you can start adding, committing and pushing changes into github by running the following commands:
  
```
git add . 
git commit -m "an informative message for you about the changes to your file(s)"
git push -u origin main
```

Anytime you will modify, cancel and/or add a new file you can check whether your directory if out of sync with github repository by running `git status`. This will prompt you a series of info telling you which file has been deleted/modified and/or added.

