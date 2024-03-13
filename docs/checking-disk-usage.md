# Checking disk usage
Below a list of commands that might be useful in case of (risking of) hitting disk quota limits. 

----------------------------------------------------------------

To check the disk quota usage within the your `HOME` and `vast/scratch` directories, do the following:

```
module load stornext/
mrquota

## This will result in something like:

--------- Vast Quota ---------
      2.974 / 100.000 TB      
    73578 / 10000000 Inodes    
--------- Home Quota ---------
       19.00 / 20.00 GB       
------------------------------
```

## Commands
* `du -hd1` -> to check the size of each file in a given folder

check