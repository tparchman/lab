# Phylogenetic analyses #

In this tutorial, we will explore how to assemble RADseq phylogenetic data sets with ``ipyrad`` and perform phylogenetic inference with quartet puzzling (``SVDquartets``/``tetrad``), maximum likelihood (``RAxML``), and Bayesian (``RevBayes``) methods. 

Email me (katuckele@gmail.com) if you have any questions or comments about this tutorial. 

## Assembly of RADseq data for phylogenetic analyses ##
``Ipyrad`` is one of the most commonly used pipelines for assembling and filtering RADseq data for phylogenetic applications. Visit the [``ipyrad`` homepage](https://ipyrad.readthedocs.io/en/latest/ "homepage") for excellent documentation and tutorials. 

### Installing ipyrad on Pronghorn ###
Installing the latest version of ``ipyrad`` is easy using conda. Detailed instructions for installation can be found [here](https://ipyrad.readthedocs.io/en/latest/3-installation.html "here"). It is recommended that you install ``ipyrad`` within the appropriate Python conda environment. At the time that this tutorial is being written, ``ipyrad`` is available for Python >=2.7 and >=3.5. 

### Create and customize the ipyrad params file ###
You can create a new ``ipyrad`` params file within Pronghorn's "base" level: 
`$ source activate py36` # activate the appropriate conda environment
`$ ipyrad -n [name of project]` # create a ipyrad params file from the command line
You can now edit your params file by providing the location of the directory containing your FASTQ files (remember to specifically target FASTQ files with "*.fastq.gz") and changing any clustering/filtering parameters. 

### Running ipyrad with sbatch script ###
Here is an example of an sbatch script that I have used to run ``ipyrad`` on Pronghorn. 

sbatch_ipyrad_test16.sh
```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --cpus-per-task 16
#SBATCH --job-name ipyrad
#SBATCH --output ipyrad_test16.txt

## activate conda environment
source activate py36

## change into the directory where your params file resides
cd /data/gpfs/assoc/parchmanlab/kuckele/jphylo

## call ipyrad on your params file
ipyrad -p params-ipyrad-test16.txt -s 1234567 -c 16 
```
You may need to change the account, partition, # CPUs, conda environment, directory, and ipyrad call as needed. 
To submit this batch script to Slurm, use the ``sbatch`` function: 
`$ sbatch sbatch_ipyrad_test16.sh`

## Quartet puzzling with tetrad (SVDquartets) ##
`Tetrad` is included in the ipyrad-analysis toolkit (see <https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-tetrad.html>), and so no additional installation is required to run `tetrad` if you have installed `ipyrad`. 

Here is an example sbatch script for running `tetrad` on Pronghorn: 

sbatch_tetrad_example.sh
```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --cpus-per-task 32
#SBATCH --job-name tetrad
#SBATCH --output tetrad_output.txt

## change into the directory where you want to print your output
cd /data/gpfs/assoc/parchmanlab/kuckele/jphylo/tetrad/

## tetrad
tetrad  -i [path to .snps.hdf5 ipyrad output file] -q 300000 -b 100 -c 32
```
This analysis will randomly sample 300,000 quartets and perform 100 bootstrap replicates. To submit the batch script to slurm, use the sbatch function: 
`$ sbatch sbatch_tetrad_example.sh`

## Maximum likelihood with RAxML ## 

### Installing RAxML on Pronghorn ###
The Exelixis Lab provides [these instructions](https://cme.h-its.org/exelixis/web/software/raxml/cluster.html "these instructions") for install on their lab's homepage. 

### Running RAxML on Pronghorn ###
Here is an example sbatch script for running RAxML on Pronghorn. 

sbatch_raxml_example.sh
```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --cpus-per-task 32
#SBATCH --job-name raxml
#SBATCH --output raxml_output.txt

## change into directory where you want RAxML output printed
cd /data/gpfs/assoc/parchmanlab/kuckele/raxml/

## call RAxML
/data/gpfs/assoc/parchmanlab/kuckele/local/bin/raxmlHPC-PTHREADS-AVX2 \
                  -T 32 \
                  -f a \
                  -m GTRGAMMA \
                  -N 100 \
                  -x 5246 \
                  -p 6554 \
                  -n raxml-jphylo \
                  -w /data/gpfs/assoc/parchmanlab/kuckele/raxml/jphyloAssemblyFull_minind4_clust85_rmchloro/ \
                  -s /data/gpfs/assoc/parchmanlab/kuckele/jphylo/jphyloAssemblyFull_minind4_rmchloro_clust85_outfiles/jphyloAssemblyFull_minind4_rmchloro_clust85_invrm.phy 
```
This above performs a rapid Bootstrap analysis (100 replicates) and searches for best-scoring ML tree. See the [RAxML manual](https://cme.h-its.org/exelixis/resource/download/NewManual.pdf "RAxML manual") for additional options. 

To submit your job to slurm, use the sbatch function: 
``$ sbatch sbatch_raxml_example.sh``

## Bayesian phylogenetic inference with RevBayes ## 
RevBayes is an interactive environment for statistical computation in phylogenetics. Multiple tutorials and example model scripts can be found at their website, <https://revbayes.github.io/>.

### Installing RevBayes on Pronghorn ###
Download the RevBayes singularity image from <https://revbayes.github.io/singularity/> into one of your home directories, e.g., 
`$ cd ~/apps/RevBayes`
`$ wget https://github.com/revbayes/revbayes/releases/download/1.1.0/RevBayes_Singularity_1.1.0.simg`
To run RevBayes, John Anderson (UNR OIT Research Computing) wrote the following sbatch script: 

rc-revbayes.sl
```
#!/bin/bash
#SBATCH  -J "revbayes"
# SBATCH -t 2:00:00
#SBATCH  --ntasks 64
#SBATCH --hint=compute_bound
#SBATCH --mem-per-cpu=3200M
#SBATCH --output=%x.%j.out              # The output file name: <job_name>.<job_id>.out
#SBATCH --error=%x.%j.err               # The error file name: <job_name>.<job_id>.err
#SBATCH --account=cpu-s5-denovo-0      # The account to charge
#SBATCH --partition=cpu-core-0   # The partition

module load openmpi/gcc/3.1.6
# Error if the input file does not exist or was not specified. Check stderr file for 
# error.
[[ -f ${1} ]] || { echo "revbayes input file does not exist" >&2; exit 1; }
# Run revbayes 
mpirun -np $SLURM_NTASKS singularity run --bind /data/gpfs/assoc --app rbmpi ~/apps/RevBayes/RevBayes_Singularity_1.1.0.simg ${1}
```
To invoke this sbatch script, use sbatch and specify your '.Rev' file, e.g., 
`sbatch rc-revbayes.sl scripts/MCMC_dating2.Rev`




