## Table of Contents
- [File System Management](#file-system-management)
  - [Checking disk quota usage](#checking-disk-quota-usage)
    - [Useful commands to manage disk quota](#useful-commands-to-manage-disk-quota)
  - [Managing data on Stornext](#managing-data-on-stornext)
    - [Archiving and retrieving files](#archiving-and-retrieving-files)
  - [VAST file system](#vast-file-system)
  - [Resources](#resources)
--------------------------------------------------

# File System Management

WEHI computer clusters (Milton) has different spaces where you can store your files. Depending on the source of data, backup requirement, and compliance requirement, you should choose different spaces. Below you will see an overview of each space on Milton (ie, `vast/projects`, `vast/scratch` and `stornext`) and how to manage your files in there.

--------------------------------------------------

## Checking disk quota usage

To check the disk quota usage, eg either within the your `HOME` and `vast/scratch/users/$USER` directories, do the following:

```
module load stornext/
mrquota
```

Which should prompt you something like:

```

--------- Vast Quota ---------
      2.974 / 100.000 TB      
    73578 / 10000000 Inodes    
--------- Home Quota ---------
       19.00 / 20.00 GB       
------------------------------
```

A `Home Quota` of ~ 20Gb will likely stop you from installing any other packages and, in case, showing you an error saying: `Disk quota exceeded`. Assuming you are not saving data in your `$HOME` directory (if so move it either to `stornext`, `vast/scratch` and/or `vast/projects`), there are many reasons why the disk quota might be full, eg:
* You are using multiple R versions each with its own packages:
* 
```
│                   
├── R
│   └── x86_64-pc-linux-gnu-library
│       ├── 4.2
|       ├── 4.3
|
```
* Some programs that you are running are saving tmp/cache data into your `$HOME` (run `ls -a $HOME` to see what's in your `$HOME` directory)
* You have created micro/mini-conda and/or mamba environments 

There are some quick checks/adjustments that you can implement to free up some space on your `$HOME` directory:

* Only have 1 R version installed 
* remove subdirectories within you `$HOME` like `.cache` which can take up lots of space. These directories will likely be re-generated when re-running the same program.
* Create your conda/mamba environments into the `vast/scratch` directory (see info on `docs/creating-conda-envs.md` for this)

### Useful commands to manage disk quota 
Below you can find a list of commands that might be useful for checking disk space and managing files on Milton:

* `du -hd1` to check the size of each file in a given folder


--------------------------------------------------

## Managing data on Stornext

Stornext is a file system that allows data to be stored and archived at the same time. Data are backed up to tape routinely which allows you to retrieve them in case of truncation. Space on stornext should be 1Tb by default, extendable to up to 5Tb but this requires justificatons from lab heads to WEHI IT support.

Regardless of the space available on stornext, you should compress each file (and eventually tar project dirs). Compressing can be done anywhere, not only on stornext, using the `pigz` module (a faster/parallel version of `gzip`) as:
```
module load pigz/ 

pigz --fast <file?> ## for faster compression 
pigz --best <file? ## slower but achieves better compression

```

### Archiving and retrieving files

Because data on stornext are backed up onto a tape-based filesystem, any given file on stornext has 2 copies: 1)a disk and 2) an archived copy. To inspect this, upon logging into Milton and navigating to the lab stornext project directory, you can run the following command:

```
module load stornext/  ## load the stornext module
snfileinfo <file> 
```

This will return you something like the below:

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
The `Location` field will tell you where the file is currently stored. From there you can see that the `snakemake.smk` file has both the `DISK AND ARCHIVE` copies. To save up space on stornext you can (and should) remove the disk copy of each file that you do not use (or plan to use) and keep the archived copy only. To do this, you need to run `snrmdiskcopy <file>`. This command will submit a job on Milton. Once the job is completed, if you run again `snfileinfo <file>` the new prompt should be:

```
-------------------------------------------------------------------------------
 File Information Report                       2024-03-22 14:09:51
 Filename:    /stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/cross-species-bulk-atac/snakefile.smk
 Stored Name: /stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/cross-species-bulk-atac/snakefile.smk
-------------------------------------------------------------------------------
      Last Modification: 15-mar-2024 16:12:44
      Owner:             unknown            Location:        ARCHIVE
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

As you can see now the `snakemake.smk` file only has the `ARCHIVED` copy, so we have freed up some disk space stornext. If you need to work with that file, you can always retrieve it back by running `snretrieve <file>`. <br/>

**NB:** Once you remove the disk copy, you can no longer access the content of that file. **however you can still remove it and/or change its name** so be careful. 

--------------------------------------------------
## VAST file system

VAST is a high-performance storage system. At WEHI, you can access `vast/scratch` and/or `vast/projects`. **No backup is available on VAST and `vast/scratch` has a 14-days deletion policy**. Each user should already have access to their own `vast/scratch` personal space, so if you type and run: `echo vast/scratch/$USER` you should see something like: `vast/scratch/vespasiani.d`. Each user should have by default 100Tb of space. <br/> 

On the other hand, to access `vast/projects` you need to complete [this form](https://support.wehi.edu.au/support/catalog/items/72), requesting the specific amount of space and for how long (up to 6 months) you'd like to have for your own project(s). Once approved, you will be able to access the project directory at `/vast/projects/project-<id>`.

An adequate and smart usage of all the 3 storage areas: 
1. `stornext`
2. `vast/scratch`
3. `vast/projects`
Should be enough to manage most projects within the lab, check the [best practices section](#best-practices) below to read how to store and work with your data.

--------------------------------------------------

## Resources
Below you can find a list of resources you can consult to learn more about filesytem management at WEHI. 

* [Data storage guidelines](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Data-how-to-store.aspx#stornext-and-vast)
* [Info on disk quotas](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Disk-quotas.aspx)
* [VAST file system](https://wehieduau.sharepoint.com/sites/rc2/SitePages/high-performance-storage.aspx)
