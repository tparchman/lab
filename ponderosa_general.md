# General information on ponderosa and associated resources

This is a group of servers for the Parchman labs private use. We currently have two independent servers, ponderosa and contorta, as well as two storage servers which are mounted with symlinks on ponderosa, where they can be accessed at `/archive` and `/mnt`.

The idea for maintaing these servers is to provide an open resource to run and or test code for specific jobs. Now that HPC has improved at UNR, pronghorn will usually be a preferred mechanism for running bigger jobs. However, using ponderosa doesnt require submiting jobs through a queueing system, so at times will be preferred. Because we all use this machine and its associated storage servers, please be sure to monitor what other are doing before starting large or time demanding jobs. The rule of thumb here is to be respectful of your lab mates, and to communicate when necessary.

## General notes on storing data

`\archive` is 38 TB of disc space in raid5 configuration. Raw sequencing data, and similar data, that needs to be stored long term should be compressed and placed at:

    /archive/parchman_lab/rawdata_to_backup/

We are backing this entire directory up to multiple other destinations, so things here should be compact. Generally all files should be gzipped or tar compressed. Please be sure to not have duplicates of these files. As we are generating more and larger sequencing runs, space gets filled up quickly, so we want to be diligent.

## General notes on running jobs

Ponderosa has 32 cores, 512 GB of RAM, and 10TB of local storage. This means read and write operations for working on ponderosa should be faster for data stored locally. For this reason, we use the directory `/working` to store data for active projects. As 10TB is not a big amount of disc space, please keep your directories within `/working` as tidy as possible:

- keep all `.fastq` files compressed whenever they are not being actively used.
- delete all `.sam` and `.bam` files when you are done processing them. These especially end up taking up enormous amounts of space.
- constantly monitor the size and content of your directories with;
    $ du -h

Every user automatically has a home directory located in `~/home/username`

## General useful commands for keeping on top of shit

To look at disc space (total, avialable, and used):

    $ df -h


