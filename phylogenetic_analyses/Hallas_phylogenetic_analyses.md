___

## RAxML
RAxML is a maxmium likelihood method used to estiamte phylogenetic relationships developed by the [Exelixis Lab](https://cme.h-its.org/exelixis/software.html), which there are variants (i.e.RAxML-Light, RAxML next generation).  For simipisicty, this will go over a basic RAxML analysis with ML tree searching and bootstrapping.  The program has many options described [HERE](https://cme.h-its.org/exelixis/resource/download/NewManual.pdf).

To use RAxML on a HPC system it is best to create a conda environment to load the program.  This elevates this issues with compling the program to allow for threading or parallelization on HPC.  *This is an issue that you will encounter if you decide to use ExaBayes, which is described below.*

```
conda create -n [name] python=[version]
```
Once an enviroment is created load the version of RAxML you wish to use.  

```
conda install -c bioconda raxml
```

I use the RAxML executable `raxmlHPC-PTHREADS`.  This allows you to use multi-core shared memory systems to speed up computational analyses.  RAxML allows you to conduct a rapid bootstrap analysis and ML tree search simultaneously, which is used with the `-f a` option.  If you choose to specify an evolutionary model that is done using the `-m` option.  There are many different models implemented in RAxML and it is best to justify the one you decide to choose.  Next we need to make sure to specifiy how many bootstraps we want to run.  This is done with the `-#` option.  I choose to use `autoMRE` instead of designating a number.  `autoMRE` determines the number of sufficient bootstrap replicates by implementing an a posteriori bootstopping analysis that uses the majority-rule consensus tree criterion to determine convergence of bootstrap analyses.  Simply, `autoMRE` option will execute a maximum of 1000 BS replicate searches, but it may, of course converge earlier.  `-x` is used as the *rapid Bootstrap Random Number Seed.*  `-p` Specify a random number seed for the parsimony inferences. This allows you to reproduce your results and will help debug the program.  `-n` name. `-o` outgroup. `-T` number of threads.  This option is only used with the `raxmlHPC-PTHREADS` executable.  *Note: if using on HPC make sure slurm resource allocation matches number of threads.  In addition, be sure to read up on the approicate number of threads.  More threads doesn't mean the anaylsis will run quicker.*

EXAMPLE
```
raxmlHPC-PTHREADS -f a -m GTRGAMMA -# autoMRE -x 12345 -p 12345 -n thamnophis_28jul.raxml.txt -s ../../ipyrad/28jul/28jul_outfiles
/28jul.phy -o NT_NT_232105 -T 16
```
Once RAxML is done, the out file *RAxML_info* file will have the number of bootstraps conducted, and the *RAxML_bipartitions* file can be used to visiualize the tree and support in figtree.
___
___

## ExaBayes
This is a bayesian analysis designed for large datasets by [Exelixis Lab](https://cme.h-its.org/exelixis/software.html), and is analogous to MrBayes or BEAST analyses.  It is fairly straight-forward but you should be come familiar with the [manual](https://cme.h-its.org/exelixis/web/software/exabayes/manual/manual.html#sec-5-2) because the program does require high volumnes of memory.  The program is also complex in it's ability for MPI. 

This program is not avaiable through conda and must be complied for the source code.  This means that you need to know which compiler and version your HPC offers.  It will have to be [downloaded](https://cme.h-its.org/exelixis/web/software/exabayes/) and transfered to HPC.

* Extract the file to working directory
```
tar -zxvf exabayes-1.5.1.tar.gz
```

* remove .tar.gz file
```
rm exabayes-1.5.1.tar.gz
```
ExaBayes needs to be complied.  It is important to know which compliers are used on the HPC.  This is because complier and openMPI are needed o be used to parallelize ExaBayes.

* Load compliers and open mpi.  Programs used on Pronghorn
```
module load openmpi/gcc/4.0.4
module load gcc/5.4.0
```

* compile from source code for MPI.
```
./configure --enable-mpi && make
```
* This program could also be compiled for threading.  The rest of the notes are in relations to the threading method.
```
./configure && make
```

There are a number of excutables that are used for ExaBayes pre and post processing.  The two main executable are `yggdrasil` and `exabayes`.  `yggdrasil` is used for threading and `exabayes` is used for MPI.  Each of these are described in detail in the manual.  You have to designated number of runs in parallele `-R` and number of chains per run `-C`.  These values also need to be stipulated in the config.nex file.  This file has all the parameter settings for the analysis.  The file is very detailed and an example file is given once ExaBayes is complied.  In this file you will also have to set how many generations you want to run ExaBayes.  

Because of the complexity of the alignment there were issues with not enough memory on BioNRES node on pronghorn.  To compensate for the lack of memory, I ran two independent runs of ExaBayes with all options used to save memory.  These options include `-M` and `-S`.   `-M` has three options: 1,2,3 and go from using more memory to less memory.  

For example, running `-R 2` and `-C 2` only finish 5% of the analysis before it was killed due to a lack of memory. Threading wasn't going to allow for enough memory to run a large dataset, and a smaller dataset was still too large.  I get a run to finish with threading I had to execute two separate analyses and combine them downstream.  I decided to use  `-R 1` and `-C 2` for 700,000 generations.  The number of generations is substationally less.  However, post evalution of this run [described below] suggested that it was roboust.

Example of bash file.
```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=jhallas@nevada.unr.edu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --mem-per-cpu=2400
#SBATCH --job-name=thamnophis_exabayes_threading_80
#SBATCH --output=output_thamnophis_exabayes_threading_80.txt

## need to load openmpi/4.0.4 before using mpi
module load openmpi/gcc/4.0.4
module load gcc/5.4.0

## change into the directory to place output files
cd /data/gpfs/assoc/denovo/jhallas/phylogenetics/thamnophis/exabayes/threading

## runs exabayes command
[full path]/yggdrasil -f [fullpath]/28jul.phy -m DNA -n myRun -s $RANDOM -c [full path]/config.nex -R 1 -C 2 -M 3 -S -T 32
```

### Post Processing

  #### exabayes output files
* ExaBayes_info.myRun
  * *Contains the same information also printed to the screen.*
* ExaBayes_topologies.myRun: 
  * *Contains all sampled topologies in nexus format*
* ExaBayes_parameters.myRun: 
  * *Contains values sampled for all non-topological, non-branch length values.  Can be viewed in Tracer*
* ExaBayes_diagnostics.myRun: 
  * *Contains chain diagnostics (e.g., acceptance ratios for all proposals or topological convergence in form of asdsf).*
    
    
Once both runs were completed I had to combine them using the following commands.  *Note: you have to make sure all files from independent runs are in the same directory*

EXAMPLE:  consense
  * creates consensus tree  

```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=2400
#SBATCH --hint=compute_bound
#SBATCH --job-name=thamnophis_consense
#SBATCH --output=output_thamnophis_consense.txt

## need to load openmpi/4.0.4 before using mpi
module load openmpi/gcc/4.0.4
module load gcc/5.4.0

## change into the directory to place output files
cd /data/gpfs/assoc/denovo/jhallas/phylogenetics/thamnophis/exabayes/threading

## runs exabayes command
/data/gpfs/assoc/denovo/exabayes-1.5.1-threading/consense -f ./ExaBayes_topologies.run* -n myCons80
```
EXAMPLE: credibleSet

*   The output file ExaBayes_credibleSet.tmp contains all sampled trees ordered by the frequency of their occurrence. You probably will find that no tree occurred more than once. 

```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=2400
#SBATCH --hint=compute_bound
#SBATCH --job-name=thamnophis_credibleSet
#SBATCH --output=output_thamnophis_credibleSet.txt

## need to load openmpi/4.0.4 before using mpi
module load openmpi/gcc/4.0.4
module load gcc/5.4.0

## change into the directory to place output files
cd /data/gpfs/assoc/denovo/jhallas/phylogenetics/thamnophis/exabayes/threading

## runs exabayes command
/data/gpfs/assoc/denovo/exabayes-1.5.1-threading/credibleSet -f ExaBayes_topologies.run* -n cred80
```
EXAMPLE: postProcParam

* Creates a  combined run table of EES and PSRF values to evaluate convergence and mixing. 

```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=2400
#SBATCH --hint=compute_bound
#SBATCH --job-name=thamnophis_postProcParam
#SBATCH --output=output_thamnophis_postProcParam.txt

## need to load openmpi/4.0.4 before using mpi
module load openmpi/gcc/4.0.4
module load gcc/5.4.0

## change into the directory to place output files
cd /data/gpfs/assoc/denovo/jhallas/phylogenetics/thamnophis/exabayes/threading

## runs exabayes command
/data/gpfs/assoc/denovo/exabayes-1.5.1-threading/postProcParam -f ExaBayes_parameters.run* -n params80
```

EXAMPLE: extractBips
* extract ESS and PSRF values for each branchlength

```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=2400
#SBATCH --hint=compute_bound
#SBATCH --job-name=thamnophis_extractBips
#SBATCH --output=output_thamnophis_extractBips.txt

## need to load openmpi/4.0.4 before using mpi
module load openmpi/gcc/4.0.4
module load gcc/5.4.0

## change into the directory to place output files
cd /data/gpfs/assoc/denovo/jhallas/phylogenetics/thamnophis/exabayes/threading

## runs exabayes command
/data/gpfs/assoc/denovo/exabayes-1.5.1-threading/extractBips -f ExaBayes_topologies.run* -n bls80
```

EXAMPLE: sdsf
*   calculates deviations of split frequencies across independent runs

```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=2400
#SBATCH --hint=compute_bound
#SBATCH --job-name=thamnophis_sdsf
#SBATCH --output=output_thamnophis_sdsf.txt

## need to load openmpi/4.0.4 before using mpi
module load openmpi/gcc/4.0.4
module load gcc/5.4.0

## change into the directory to place output files
cd /data/gpfs/assoc/denovo/jhallas/phylogenetics/thamnophis/exabayes/threading

## runs exabayes command
/data/gpfs/assoc/denovo/exabayes-1.5.1-threading/sdsf -f ExaBayes_topologies.run*
```
___
___

## Multiple Species Coalescent model in tetrad

This similar to SVDquartet.  Multiple Species Coalescent model is reviewed in Edwards et al., 2016; Liu et al., 2019.  Unlike other methods that take a two-step approach by first generating many individual gene trees and then summarizing the resulting gene trees to infer the species tree (e.g. ASTRAL (Zhang et al., 2018), MP-EST (Liu et al., 2010), or NJst (Liu and Yu, 2011)), tetrad utilizes unlinked SNPs from aligned sequence data to generate quartets under a coalescent model then uses the wQMC algorithm to join the resulting quartets into a species tree. 

Manual for tetrad located [here](https://github.com/eaton-lab/tetrad).  `-i` input file. The number of cores `-c`.  number of bootstraps `-b`.  Number of quartets`-q`. A very large number will evaluate all quartets.  Output file will report this.   `-n` is name for output file

tetrad uses .snps.hdf5 generated from ipyrad 

load tetrad through conda
```
conda install -c eaton-lab tetrad 
```


EXAMPLE
```
#!/usr/bin/env bash
#SBATCH -A cpu-s1-bionres-0
#SBATCH -p cpu-s1-bionres-0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=jhallas@nevada.unr.edu
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --hint=compute_bound
#SBATCH --mem-per-cpu=2400
#SBATCH --job-name thamnophis_tetrad_3
#SBATCH --output output_thamnophis_26Oct_tetrad_3.txt

source activate py36

## change into the directory where your params file resides
cd /data/gpfs/assoc/denovo/jhallas/phylogenetics/thamnophis/tetrad

## call ipyrad on your params file
tetrad -i 28jul.snps.hdf5 -c 32 -b 50 -q 1000000 -n [name]
```



