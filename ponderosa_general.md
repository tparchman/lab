# General information on ponderosa and associated resources

This is a group of servers for the Parchman labs private use. We currently have three independent servers (ponderosa, contorta, and wallace) , as well as two storage servers which are mounted with symlinks on ponderosa, where they can be accessed at `/backups` and `/mnt`.

The idea for maintaing these servers is to provide an open resource to run and or test code for specific jobs. Now that HPC has improved at UNR, pronghorn will usually be a preferred mechanism for running bigger jobs. However, using ponderosa doesn't require submiting jobs through a queueing system, so at times will be preferred. Because we all use this machine and its associated storage servers, please be sure to monitor what other are doing before starting large or time demanding jobs. The rule of thumb here is to be respectful of your lab mates, and to communicate when necessary.

## General notes on storing data IGNORE UNTIL UPDATE

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


# Two major rules:

 ## 1. **keep all files compressed** when you are not actively using them

 ## 2. **rm all files that are not absolutely necessary**


## General useful commands for keeping on top of things
See thorough cheat sheet at the bottom for other stuff.


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

### Example user add, for Abby Miller

Use `sudo` to activate account, set working directory

    $ sudo useradd -m -s /bin/bash -c "Angie Lenard, Parchman Group" -G users,working alenard

this adds a new user, alenard:

- `m` creates home directory and copies files from /etc/skel
- `s` /bin/bash: makes bash the default shell
- `c` "Angie Lenard, Parchman Group" adds comment to /etc/passwd file
- `G` users, working adds user to secondary group working.

Set passwd:
 
    $ sudo passwd alenard

 
Set the passwd to G00gle_it (temporary)

Age password so user will have to change the first time they login

    $ sudo chage -m 10 alenard

Then to login: 

    $ ssh alenard@ponderosa.biology.unr.edu
    password: G00gle_it (temporary; those are zeros not ones.)
    
    change password during first login using:

    $ passwd <newpassword>

# Linux Command Cheat Sheet

This cheat sheet summarizes the commands, grouped by topic.

---

## Shell Basics & Navigation
| Command | Purpose | Example |
|---------|---------|---------|
| `pwd` | Print current working directory | `pwd` |
| `ls` | List files in directory | `ls -F` |
| `ll` | Long format list (alias for `ls -laF`) | `ll` |
| `cd` | Change directory | `cd ~/Documents` |

---

## Shell Customization (Zsh/Oh My Zsh)
| Command | Purpose | Example |
|---------|---------|---------|
| `brew install zsh` | Install Zsh on macOS | `brew install zsh` |
| `sudo apt install zsh` | Install Zsh on Ubuntu | `sudo apt install zsh -y` |
| `chsh -s $(which zsh)` | Make Zsh the default shell | `chsh -s $(which zsh)` |
| `alias` | Create shortcut for a command | `alias ll="ls -laF"` |
| `set -o noclobber` | Prevent overwriting with `>` | `set -o noclobber` |

---

## Package Managers
| Command | Purpose | Example |
|---------|---------|---------|
| `sudo apt update` | Update package list (Ubuntu) | `sudo apt update` |
| `sudo apt upgrade` | Upgrade installed packages | `sudo apt upgrade` |
| `sudo apt install <pkg>` | Install a package | `sudo apt install curl` |
| `brew update` | Update Homebrew package list (macOS) | `brew update` |
| `brew upgrade` | Upgrade installed brew packages | `brew upgrade` |
| `brew install <pkg>` | Install a package | `brew install wget` |

---

## File Management
| Command | Purpose | Example |
|---------|---------|---------|
| `cp` | Copy files | `cp file1.txt file2.txt` |
| `mv` | Move/rename files | `mv old.txt new.txt` |
| `rm` | Remove files (interactive if aliased) | `rm file.txt` |
| `rm -rf` | **Dangerous**: force remove recursively | `rm -rf data/` |
| `rmdir` | Remove empty directory | `rmdir olddir` |
| `split` | Split file into chunks | `split -l 1000 big.fastq chunk_` |
| `cat` | Concatenate/display files | `cat chunk_* > combined.fastq` |
| `head` | Show first *n* lines | `head -n 4000 sample.fastq` |
| `tail` | Show last *n* lines | `tail -n 100 sample.fastq` |

---

## Text Processing & Pipes
| Command | Purpose | Example |
|---------|---------|---------|
| `grep` | Search by pattern | `grep "chrIII" yeast.gff` |
| `grep ^@ file` | Match lines starting with `@` | `grep ^@ sample.fastq` |
| `tr` | Translate/replace characters | `tr 'T' 'U' < seqs.txt` |
| `sort` | Sort lines | `sort features.txt` |
| `uniq` | Remove duplicates (after `sort`) | `sort features.txt \| uniq` |
| `cut` | Extract columns from tab-delimited file | `cut -f3 yeast.gff` |
| `wc -l` | Count lines | `grep ^@ sample.fastq \| wc -l` |
| `|` (pipe) | Send output to another command | `ls \| wc -l` |
| `>` | Redirect output (overwrite) | `echo "hi" > out.txt` |
| `>>` | Redirect output (append) | `echo "again" >> out.txt` |

---

## Process Control & Job Management
| Command | Purpose | Example |
|---------|---------|---------|
| `top` | Monitor processes in real time | `top` |
| `htop` | Interactive process monitor | `htop` |
| `ps` | List your processes | `ps` |
| `ps aux` | List all processes with details | `ps aux` |
| `ps aux \| grep <name>` | Find process by name | `ps aux \| grep python` |
| `pgrep -a <name>` | Find PID(s) by process name | `pgrep -a ping` |
| `kill <PID>` | Kill process by ID | `kill 12345` |
| `pkill <name>` | Kill process by name | `pkill yes` |
| `ctrl c` | Stop a running foreground job | *(keyboard shortcut)* |
| `ctrl z` | Suspend job | *(keyboard shortcut)* |
| `bg` | Resume suspended job in background | `bg` |
| `<cmd> &` | Run job in background | `yes > /dev/null &` |
| `nohup <cmd> &` | Run job immune to hangups | `nohup ping google.com &` |

---

## Number Generators
| Command | Purpose | Example |
|---------|---------|---------|
| `jot -r N` | Generate N random numbers (macOS) | `jot -r 100` |
| `seq` | Generate sequences of numbers | `seq 1 10` |
| `shuf` | Generate random numbers (Linux) | `shuf -i 1-100 -n 10` |

---

## Permissions
| Command | Purpose | Example |
|---------|---------|---------|
| `ls -l` | List files with permissions | `ls -l` |
| `chmod a+x file.sh` | Make executable | `chmod a+x script.sh` |
| `chmod a+w file` | Add write permission | `chmod a+w file.txt` |

---

## Remote Access & Transfers
| Command | Purpose | Example |
|---------|---------|---------|
| `ssh user@server` | Log into remote server | `ssh user@hpc.edu` |
| `sftp user@server` | Transfer files interactively | `sftp user@hpc.edu` |
| `rsync -av src/ dest/` | Copy/sync directories | `rsync -av data/ backup/` |
| `rsync -av --delete src/ dest/` | Sync exactly (delete removed files) | `rsync -av --delete data/ backup/` |
| `rsync -av user@server:/src/ dest/` | Copy from remote server | `rsync -av user@hpc.edu:/scratch data/` |
| `curl <url> -o file` | Download from URL | `curl "https://.../seq.fasta" -o seq.fasta` |
| `wget <url>` | Download from URL | `wget https://.../yeast.gff` |

---

## Shell Scripting
| Command | Purpose | Example |
|---------|---------|---------|
| `#!/bin/bash` | Shebang line for bash scripts | *(first line in script)* |
| `#!/bin/zsh` | Shebang line for zsh scripts | *(first line in script)* |
| `echo` | Print text to screen | `echo "Hello world"` |
| `bash script.sh` | Run script with bash | `bash script.sh` |
| `chmod +x script.sh` | Make script executable | `chmod +x script.sh` |
| `./script.sh` | Run executable script | `./script.sh` |

---

## Safety Aliases (from our `.zshrc`)
```zsh
alias python='python3'
alias ll='ls -laF'
alias ls='ls -F'
alias rm='rm -i'
alias mv='mv -i'
alias cp='cp -i'
