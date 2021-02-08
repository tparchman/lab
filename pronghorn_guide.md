## General notes on pronghorn usage

### Where to work on pronghorn

For most of parchman lab work, things should be stored and worked on from:

    /data/gpfs/assoc/parchmanlab/

Files can be easily moved to or from here using scp or rsync.

    $ scp tparchman@pronghorn.rc.unr.edu:/data/gpfs/assoc/parchmanlab/parchman/*.pl .

    $ rsync -av muricata_fastqs tparchman@pronghorn.rc.unr.edu:/data/gpfs/assoc/parchmanlab/parchman/

### Containers or environments for managing software

Containers or environments must be used on pronghorn to install and run user specific software. Using such tools means that you can control your own portable environment, and that individual users do not have to install software more generally on the system. In the parchman lab we have been using `anaconda` environments for this. Many others are using `singularity` containers.

%difference between miniconda and anaconda

### Install `anaconda` on pronghorn

Using `anaconda` will allow you to install all of the programs you need on pronghorn. Build within user-specific directory, e.g.: 
    
    /data/gpfs/assoc/parchmanlab/parchman/
	
    $ wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
    
    $ bash Anaconda3-2020.11-Linux-x86_64.sh 

It will make the most sense to put anaconda3 in the directory that serves as your main working directory.

    $ rm Anaconda3-2019.10-Linux-x86_64.sh

To activate conda either logout/log back in or:

    $ source .bashrc

### Create your main working environment and add some channels

    $ conda create -n py38 python=3.8

%describe `channels`

    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda config --add channels defaults
    $ conda config --add channels r

NOTE: may need to modify .condarc in order to move conda-forge and bioconda to top of list

## activating and deactivating environment

%need a few sentences describing [base], how to deactivate out of base, how to activate back into base, the difference between base and your environment.

%warnings about installing in base, and why you dont want to do that.

To activate your environment, use

    $ conda activate py38

To deactivate an active environment, use

    $ conda deactivate

Note, if you deactivate after activating your named environment, you can deactivate out of base using the same command.

## conda software installs

Paragraph about considerations before choosing which version of what to install.

To figure out which channel and which version of a given software you want to install, it is generally best to google what you are looking for to see what is available and what channel you want to use.

Below will install the most current version of bwa available from `bioconda`

    $ conda install -c bioconda bwa
    $ conda install -c bioconda samtools
    $ conda install -c bioconda bcftools
    $ conda install -c bioconda vcftools

## conda software uninstalls

## removing anaconda environment

    $ conda install anaconda-clean
    $ anaconda-clean --yes

### Building additional environments to keep dependencies clean

 You will likely want to build different conda environments for different software installs to keep dependencies clean. `Ipyrad` is an example of where we have built separate environments.