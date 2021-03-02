# General information on ponderosa and associated resources

This is a group of servers for the Parchman labs private use. We currently have two independent servers, ponderosa and contorta, as well as two storage servers which are mounted with symlinks on ponderosa, where they can be accessed at `/archive` and `/mnt`.

The idea for maintaing these servers is to provide an open resource to run and or test code for specific jobs. Now that HPC has improved at UNR, pronghorn will usually be a preferred mechanism for running bigger jobs. However, using ponderosa doesnt require submiting jobs through a queueing system, so at times will be preferred. Because we all use this machine and its associated storage servers, please be sure to monitor what other are doing before starting large or time demanding jobs. The rule of thumb here is to be respectful of your lab mates, and to communicate when necessary.

## General notes on storing data

`\archive` is 38 TB of disc space in raid5 configuration. Raw sequencing data, and similar data, that needs to be stored long term should be compressed and placed at:

    /archive/parchman_lab/rawdata_to_backup/

We are backing this entire directory up to multiple other destinations, so things here should be compact. Generally all files should be gzipped or tar compressed. Please be sure to not have duplicates of these files. As we are generating more and larger sequencing runs, space gets filled up quickly, so we want to be diligent.

## General notes on running jobs

