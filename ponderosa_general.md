# General information on ponderosa and associated resources

This is a group of servers for the Parchman labs private use. We currently have two independent servers, ponderosa and contorta, as well as two storage servers which are mounted with symlinks on ponderosa, where they can be accessed at `/archive` and `/mnt`.

The idea for maintaing these servers is to provide an open resource to run and or test code for specific jobs. Now that HPC has improved at UNR, pronghorn will usually be a preferred mechanism for running bigger jobs. However, using ponderosa doesn't require submiting jobs through a queueing system, so at times will be preferred. Because we all use this machine and its associated storage servers, please be sure to monitor what other are doing before starting large or time demanding jobs. The rule of thumb here is to be respectful of your lab mates, and to communicate when necessary.

## General notes on storing data

`/archive` is 38 TB of disc space in raid5 configuration. Raw sequencing data, and similar data, that needs to be stored long term should be compressed and placed at:

    /archive/parchman_lab/rawdata_to_backup/

We are backing this entire directory up to multiple other destinations, so things here should be compact. Generally all files should be gzipped or tar compressed. Please be sure to not have duplicates of these files. As we are generating more and larger sequencing runs, space gets filled up quickly, so we want to be diligent.

## General notes on running jobs

Ponderosa has 32 cores, 512 GB of RAM, and 10TB of local storage. This means read and write operations for working on ponderosa should be faster for data stored locally. For this reason, we use the directory `/working` to store data for active projects. As 10TB is not a big amount of disc space, please keep your directories within `/working` as tidy as possible:

- keep all `.fastq` files compressed whenever they are not being actively used.
- delete all `.sam` and `.bam` files when you are done processing them. These especially end up taking up enormous amounts of space.
- constantly monitor the size and content of your directories with:
```
    $ du -h
```
Every user automatically has a home directory located in `~/home/username`

## Software installs and modules

Software installed and compiled at the system level by system administrator (Mike Zierten, mzierten@unr.edu) or Parchman is maintained in a module library. If you need something installed that is not already at the system level, contact Tom and Mike Zierten, and forward a link to the specific software you need installed. Mike gets pulled on by a lot of people in a lot of different directions, so keep that in mind if you ever need his help. He is an extremely nice guy, but please be organized to save him time, and make sure he knows how much we appreciate him.

To view software modules available:

    $ module avail

To load a module so that you can use the software:

    $ module load name_of_module

For example, to load bwa/0.7.8

    $ module load bwa/0.7.8
    
## Emergencies?

Anytime the system appears to be down, or you can not log in, immediately contact Tom or, if necessary, Mike Zierten. Mike is the College of Science system administrator for research computing, and he helps us keep things stable.

## General useful commands for keeping on top of shit

To look at disc space (total, avialable, and used):

    $ df -h

Monitor disc usage

    $ du -h 

To check the usage for each user within a certain directory

    $ du -sh *

Simple file compression with `gzip`

    $ gzip data.fastq

Simple file decompression with `gzip`

    $ gunzip data.fastq.gz

Example of using tar to compress directory:

    $ tar -zcvf directory.tgz directory/

Example of using tar to decompress:

    $ tar -zxvf directory.tgz

### More specific, but commonly useful commands

running jobs in background without interupt
    
    $ nohup <COMMANDS> &>/dev/null &
    
changing directory permissions for one user at a time:

    $ setfacl -m u:jpjahner:rwx test/

Renaming a batch of files

    $ rename 's/out_//g' *

Directory trees with full info

    $ tree -s -h directory_of_choice/

Killing all jobs by user (only for those with sudo)
    
    $ sudo pkill -u faske

## Setting up global protect VPN for off campus ssh access

This is a new requirement as of December of 2022. Before ssh connection to ponderosa, VPN through global protect needs to be activated. All ssh connections from off campus need to occur **after** VPN connection is established.

**VPN connection prior to ssh is only necessary if you are connecting from offcampus.** If you are on campus, no need to worry about this.

### Globalprotect VPN client needs to installed

download for mac or win at https://border.unr.edu/global-protect/getsoftwarepage.esp

### Access to global connect must be requested from OIT

The below link will take you where you need to go.

https://unr.teamdynamix.com/TDClient/2684/Portal/KB/ArticleDet?ID=117539
 
The link to requesting VPN access is the line that says “granted access to the VPN”. All that is required for access permission is a UNR netID.

### Using Globalprotect

Go to connect portal, and type in the below address:

    vpn.unr.edu

Enter Netid, then password. Multifactor authentication (MFA) will then occur via text or email. Once you are connected, you can leave it on, and `ssh` at will. 

You only should only need to enter your netID and password once, which makes using this app exceptionally easy. You can leave connected, or just reconnect before using ssh. This should not have any consequences for connection speed or internet usage.



## Setting up ssh so that you dont have to use a password

When you `ssh` to a remote server, e.g., ssh tparchman@pronghorn.rc.unr.edu, you are prompted for a password. Because you are doing this many times a day, this is annoying. It is also not necessary. By carefully creating, `id_rsa` and `id_rsa.pub` in `.ssh`, and then appending your `id_rsa.pub` key to the `authorized_keys` file in `.ssh` on the remote server, you set up these keys so that the server recognizes login attempts from your local computer. id_rsa contains an id that identifies your machine. `id_rsa.pub` contains the public id to be recognized for remote servers.

1. Go to `.ssh` on your private computer (or whichever computer will be the login source). You should find `id_rsa` and `id_rsa.pub`. If they are there, have a look so that you understand the information they store.

	If these files do not exist, you need to make them.

	    $ ssh-keygen -t rsa

	Generating public/private rsa key pair. 

	IMPORTANTLY, you dont need to create a passphrase, and I suggest not doing so. When prompted for passphrase just hit ENTER for blank passphrase.

2. On the remote server (ponderosa, pronghorn, or wherever else), go to .ssh in your home directory. You should find a file, `authorized_keys`, which contains keys for whichever computers you choose. Have a look so that you understand the information it contains. IF it doesn't exist, create a blank file with this name.


	    $ less authorized_keys
	
    or
    
	    $ touch authorized_keys

3. Copy the `id_rsa.pub` file from you local computer to `.ssh` on the remote server (use care if id_rsa.pub already exists on the server, in which case I'd suggest copying `id_rsa.pub` on your local computer into a different file name before moving).

	    $ scp id_rsa.pub username@pronghorn.rc.unr.edu:~/.ssh/

4. Append the contents of `id_rsa.pub` (from your local computer) to `authorized_keys` (on the remote server).

	    $ cat id_rsa.pub >> authorized_keys

	or just copy the line of interest in id_rsa.pub and copy it into authorized keys using nano, emacs, or vim.

## For sudo only: setting up user accounts

### Example user add, for Austin Baldridge

Use `sudo` to activate account, set working directory

    $ sudo useradd -m -s /bin/bash -c "Ellen McMullen, Parchman Group" -G users,working emcmullen

this adds a new user, emcmullen:

- `m` creates home directory and copies files from /etc/skel
- `s` /bin/bash: makes bash the default shell
- `c` "Ellen Mcmullen, Parchman Group" adds comment to /etc/passwd file
- `G` users, working adds user to secondary group working.

Set passwd:
 
    $ sudo passwd emcmullen
 
Set her passwd to G00gle_it (temporary)

Age passwork so user will have to change the first time they login

    $ sudo chage -m 10 emcmullen

Then to login: 

    $ ssh emcmullen@ponderosa.biology.unr.edu
    password: G00gle_it (temporary; those are zeros not ones.)
    
    change during first login using:

    $ passwd newpassword

