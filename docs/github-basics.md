# Github basics
----------------------------------------------------------------
## Table of Contents
- [Github basics](#github-basics)
  - [Table of Contents](#table-of-contents)
  - [Configuration](#configuration)
    - [Setting up SSH key](#setting-up-ssh-key)
  - [Basic git usage](#basic-git-usage)

----------------------------------------------------------------

## Configuration

### Setting up SSH key
The below information on how to set up github push with ssh keys come from [this post](https://gist.github.com/xirixiz/b6b0c6f4917ce17a90e00f9b60566278) and [this github page](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)

Create a new ssh key to be linked with the github email address by running this command either on Windows or macOS: `ssh-keygen -t rsa -b 4096 -C "genomicskinglab@gmail.com"` and press `Enter` when the terminal will prompt the following:
```
Enter file in which to save the key (~/.ssh/id_rsa):
Enter passphrase (empty for no passphrase):  
Enter same passphrase again:
```
Copy the newly generated key on your keyboad by running `pbcopy < ~/.ssh/id_rsa.pub` on macOS or `cat ~/.ssh/id_rsa.pub | clip` on Windows. **NB:** this is just like running `CTRL + C` on a selected text. Then go into the lab's github.com account > Settings > SSH and GPG keys, click onto the `New SSH key` button, copy the key in there (ie `CTRL + V`) and give it a meaningful title (eg `<your-name>-ssh-authentication-key`). <br/>

Next, go back to your terminal and run `ssh -T git@github.com` to check whether the ssh key successfully works and you can be authenticated into the lab's github account. The terminal should prompt you something like this:
```
Hi Genomics-KingLab! You've successfully authenticated, but GitHub does not provide shell access.
```

Once you successfully connected your local machine to github you can run `git add` , `git commit -m 'your message'` and `git push` to push your changes onto github. Finally, add the ssh key to the agent by running `ssh-add ~/.ssh/id_rsa`. **NB:** this might not work on Windows (not sure why).


## Basic git usage

Read [this post](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Basic-Git-Usage.aspx) and [this repo](https://github.com/WEHI-ResearchComputing/GitIntro/blob/master/Lessons/GitIntro_1.md) from WEHI IT support to learn how to use the basic git functionalities. 