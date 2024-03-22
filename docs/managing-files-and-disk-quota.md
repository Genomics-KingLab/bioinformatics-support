# Managing files on Stornext

Check [this documentation](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Data-how-to-store.aspx#stornext-and-vast) to understand how filesystem management works on Milton.

## Archiving and retrieving files

Files saved on stornext are automatically backed up on a tape-based filesystem. Each file on stornext automatically has 2 copies: 1) disk and 2) archive. To manage files and disk quota on Milton you need to load the stornext module as `module load stornext/`. <br/>

To check the status of each file within stornext you need to run `snfileinfo <file>`. This will return you something like the below:
```
-------------------------------------------------------------------------------
 File Information Report                       2024-03-22 13:59:31
 Filename:    /stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/cross-species-bulk-atac/snakefile.smk
 Stored Name: /stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/cross-species-bulk-atac/snakefile.smk
-------------------------------------------------------------------------------
      Last Modification: 15-mar-2024 16:12:44
      Owner:             unknown            Location:        DISK AND ARCHIVE
      Group:             unknown            Existing Copies: 2
      Access:            644                Target Copies:   2
                                            Expired Copies:  0
      Target Stub:       0 (KB)             Existing Stub:   n/a (KB)
      File size:         3,277              Store:           MINTIME
      Affinity:          n/a                Reloc:           MINTIME
      Class:             general_policy     Trunc:           MINTIME
      Alt Store Copy:    Disabled           Clean DB Info:   NO
      Media:      PC0552(1)  UB0600(2)  
      Checksum:   N
      Encryption: N
      Object Ids: N
```

To save up space on stornext you can (and should) remove the disk copy of each file that you do not use and keep the archived copy only. To do this, you need to run `snrmdiskcopy <file>`. Once the job is completed, if you run again `snfileinfo <file>` the new prompt should be:

```
-------------------------------------------------------------------------------
 File Information Report                       2024-03-22 14:09:51
 Filename:    /stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/cross-species-bulk-atac/snakefile.smk
 Stored Name: /stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/cross-species-bulk-atac/snakefile.smk
-------------------------------------------------------------------------------
      Last Modification: 15-mar-2024 16:12:44
      Owner:             unknown            **Location:        ARCHIVE**
      Group:             unknown            Existing Copies: 2
      Access:            644                Target Copies:   2
                                            Expired Copies:  0
      Target Stub:       0 (KB)             Existing Stub:   0 (KB)
      File size:         3,277              Store:           MINTIME
      Affinity:          n/a                Reloc:           MINTIME
      Class:             general_policy     Trunc:           MINTIME
      Alt Store Copy:    Disabled           Clean DB Info:   NO
      Media:      PC0552(1)  UB0600(2)  
      Checksum:   N
      Encryption: N
      Object Ids: N

```
Once you removed the disk copy, you can no longer access the content of that file. **NB: you can still remove it and/or change its name**. In case you need to work with that file, you can retrieve it back by running `snretrieve <file>`.


## Checking disk usage

To check the disk quota usage within the your `HOME` and `vast/scratch/users/$USER` directories, do the following:

```
module load stornext/
mrquota

## The command should prompt you something like:

--------- Vast Quota ---------
      2.974 / 100.000 TB      
    73578 / 10000000 Inodes    
--------- Home Quota ---------
       19.00 / 20.00 GB       
------------------------------
```

## Useful commands to manage the disk quota 
Below you can find a list of commands that might be useful for checking disk space and managing files on Milton
* `du -hd1` to check the size of each file in a given folder

