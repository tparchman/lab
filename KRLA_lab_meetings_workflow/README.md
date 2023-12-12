# Genotyping By Sequencing (GBS) - An Example Study

This document presents our workflow and rationale for genotype inference from high throughput sequencing of reduced representation libraries (GBS, RADseq, ddRADseq, etc). Several canned software packages or computational workflows exist for handling this type of data. These methods rely on set thresholds for sequencing coverage depth per locus to call hard genotypes. The biggest cost of using these methods is throwing away much, if not most, of the data.

The methods described in this document are meant to be thorough and allow for user control at each step of the bioinformatic process. However, 'simpler' alternatives exist in the form of canned software packages and other streamlined computational workflows. These methods all rely on set thresholds for sequencing coverage depth per locus to call hard genotypes. This is highly problematic. Biggest cost of using these types of methods is throwing away much if not most of the data. Examples of such software/workflows include:
* [Stacks](http://catchenlab.life.illinois.edu/stacks/)
* [Ddocent](http://www.ddocent.com//)
* [ipyrad](https://ipyrad.readthedocs.io/en/master/)

See [Nielsen et al. 2011](./papers/Nielsen_etal_2011.pdf) and [Buerkle and Gompert 2013](./papers/Buerkle_Gompert_2013.pdf) for articulate thoughts about this.

## GBS Workflow Table of Contents

1) [INITIAL SEQUENCE PROCESSING](#1-initial-sequence-processing)\
    a. [Contaminant cleaning using tapioca](#1a-cleaning-contaminants)\
    b. [Parsing barcodes](#1b-barcode-parsing)\
    c. [Splitting fastqs](#1c-splitting-fastqs)
2) [DENOVO REFERENCE ASSEMBLY](#2-denovo-assembly-to-generate-a-consensus-reference-for-mapping-reads-prior-to-genotyping)\
    a. [Directory & file prep](#2a-directory--file-prep)\
    b. [Generating unique sequence files](#2b-generate-unique-sequence-files-for-each-individual)\
    c. [Sequence subsetting for alignment](#2c-subset-sequences-for-contig-alignment-and-assembly)
3) [READ MAPPING](#3-mapping-reads-from-all-individuals-to-reference-using-bwa)\
    a. [Directory & file prep]()
4) [CALLING VARIANTS]()
5) [FILTERING]()
6) [GENOTYPE PROBABILITIES]()

**... add more as we work**

### To do:

* Add explanation of selectContigs.sh parts
* Add estimated run times for each step
* Add info about how fastq / assmebly folders are used? 
* Add full paths for scripts used

## Species information (*Krascheninnikovia lanata*) 
<img src="images/KRLAplant.jpg" width="310"> &emsp; &emsp;
<img src="images/KRLAmap.png" width="401">

&emsp; *Krascheninnikovia lanata* (winterfat) is a perennial shrub with a broad north/south distributional range spanning western Canada, U.S. and Mexico. The species is exclusive to North America and its current range is likely the result of southward expansion following two distinct migration events from eastern Mongolian lineages ~ 1.8 - 0.5 Mya.\
&emsp; The species is a halophyte (salt-tolerant) and is one of the only species outside of the *Atriplex* complex to co-dominate the salt desert shrublands of the Great Basin. It is a highly nutritious source of forage which is notable given that it is prone to replacement by the toxic exotic *Halogeton glomeratus* within disturbed habitats. The common name 'winterfat' is indicative of the persistence of green leaves throughout the winter season and late fall phenology.\
&emsp; Population sampling was a combined effort throughout 2021 - 2022 with Cathy Silliman doing collections for most of the populations in the west, central, and north Great Basin and Seth Romero gathering collections from the eastern Great Basin and Mojave. Anecdotally, many of the individuals in the north, central and eastern Great Basin were smaller in stature but part of broad near-monocultures that created consistent cover across broad areas. By contrast, many of the populations in the Mojave and southwest Great Basin were large individuals that occured in small islands or as sub-dominants with only a handful of individuals living in close proximity. The populations **DT** and **CL**, in particular, had nearly every individual sampled that could found at those locations.

## File structure

* Create a project folder (wherever you want it) called KRLA with `mkdir KRLA`
* Subfolders will be created throughout the analysis

## Sample collection

* Individuals were sampled across their natural distribution
* Information about sampling can be found at [insert doc path here]

## DNA extraction

* Performed by AG Biotech in 12/2022
* DNA information with plate maps and IDs can be found at [insert doc path here]

## Library preparation

* Performed in Parchman Lab in 12/2022
* [add upstream path]KRLA\_RFseq\_mastermixcockatils.xlsx contains information about reagents used in library prep
* R/L and PCR for plates 1-6

## Sequencing

* 1 lane of S2 chemistry NovaSeq data at UTGSAF in 03/2023
* sequence results in [add upstream path]KRLA\_S1\_L001\_R1\_001.fastq.gz

## Data cleaning

### Cleaning contaminants

* **Goal:** Remove reads that are not from the target organims
* Performed 06/09/23

1. Create a clean\_data folder inside KRLA folder:
	
	```sh
	mkdir clean_data
	```
	
2. Copy the sequencing results to your clean\_data folder:

	```sh
	cd clean_data
	cp /archive/parchman_lab/rawdata_to_backup/KRLA/KRLA_S1_L001_R1_001.fastq.gz .
	```
	
3. Decompress the file with:

	```sh
	gunzip KRLA_S1_L001_R1_001.fastq.gz
	```

4. Count the reads **before** cleaning with:

	```sh
	$ nohup grep -c "^@" KRLA_S1_L001_R1_001.fastq > KRLA_number_of_rawreads.txt &
	```
	
	* `nohup` (no hang up) allows the command to keep running, even if you close the terminal, lose connection with the server, or log out.
	* `&` allows the command to run in the background.
	* `grep -c` will count the number of occurances of the pattern ("^@") in the file (KRLA_S1_L001_R1_001.fastq).
	* the result file (KRLA_number_of_rawreads.txt) will contain the initial number of reads 
	* **Reads before cleaning:** n

5. Run the cleaning script (note: the script will need to be edited to change the paths to the appropriate files):

	```sh
	cp path/to/cleaning/script.sh .
	nano cleaning_bash_KRLA.sh
	conda deactivate
	module load fqutils/0.4.1
	module load bowtie2/2.2.5
	bash cleaning_bash_KRLA.sh &
	```
	
	* `module load` is needed to allow access to the modules used in the script
	* fqutils and bowtie2 are used in the cleaning script
	* cleaning\_bash\_KRLA.sh uses the tapioca pipeline's contamination analysis script to remove artifacts from:
		* Illumina oligos
		* PhiX
		* *E. coli*

6. Check that KRLA.clean.fastq has been created (from cleaning script), then remove the raw data file (from KRLA/data - a backup is stored at /archive/parchman\_lab/rawdata\_to\_backup/KRLA) :

	```sh
	ls
	rm -rf KRLA_S1_L001_R1_001.fastq &
	```
	
7. Count the reads **after** cleaning with:

	```sh
	nohup grep -c "^@" KRLA.clean.fastq > KRLA_clean_reads.txt &
	```
	
	* reads after cleaning are typically ~75% of reads before cleaning
	* **Reads after cleaning:** n


### Barcode parsing

* **Goal:** Remove barcode and adapter sequences from reads, place this information (what individual the reads came from) in the read ID line

1. Run parsing script:

	```sh
	cp /working/parchman/KRLA/parse_barcodes769.pl .
	cp /working/parchman/KRLA/KRLA_barcode_key.csv .
	nohup perl parse_barcodes769.pl KRLA_barcode_key.csv KRLA.clean.fastq A00 &>/dev/null &
	```
	
	* KRLA\_barcode\_key.csv provides the barcodes for each individual
	* KRLA.clean.fastq are the clean reads from contaminant cleaning
	* A00 is the sequencer ID (first 3 characters after @ in the fastq identifier)
	* &>/dev/null prevents display of stdout

2. View the parsing results (created by the parsing script):

	```sh
	less parsereport_KRLA.clean.fastq
	```
	
	* results compare good and bad mids (molecular IDs), bad mids are usually <5% of total mids
	* the other removed results are typically insignificant
	* **Parse report:**
		* Good mids count: 1571963061
		* Bad mids count: 73508410
		* Number of seqs with potential MSE adapter in seq: 321195
		* Seqs that were too short after removing MSE and beyond: 428

### Splitting fastqs

* **Goal:** Sort the reads that are in the one .fastq file into multiple .fastq files, one for each individual

1. Create an IDs file:

	```sh
	cut -f 3 -d "," KRLA_barcode_key.csv | grep "_" > KRLA_ids_noheader.txt
	```
	
	* `-f 3` option looks at the 3rd field (column)
	* `-d ","` option specifies comma as delimiter
	* `cut` extracts that column and pipes it to `grep`
	* `grep` ensures the ID column extracted contains an underscore (this removes the header line / label)

2. Split the (cleaned and parsed) fastq files by individual ID:

	```sh
	cp /working/parchman/KRLA/splitfastqs/splitFastq_universal_regex.pl .
	nohup perl splitFastq_universal_regex.pl KRLA_ids_noheader.txt parsed_KRLA.clean.fastq &>/dev/null &
	```
	
	* the perl script creates a separate fastq file for each individual's reads

3. Compress all of the resultant fastq files:

	```sh
	nohup gzip *fastq &>/dev/null &
	```

## Denovo assembly

Denovo assemblies uses your reads to generate a consensus artificial reference to map your reads against. This is the only option in some systems where there is no reference genome available. Sometimes it is still prefered, depending on the genome divergence between your species or population and that of the reference genome. If you are using a reference genome, download the reference and skip this step.

There are many tools available for denovo assembly. We will be using [CD-HIT](https://www.bioinformatics.org/cd-hit/) in a workflow from [dDocent](http://www.ddocent.com/) (Puritz et al. 2014). See [LaCava et al. 2020](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13108) for a comparison of tools.

This workflow begins with the gziped .fastq files in KRLA/clean\_data.

### Prepare directories and files

1. Create subfolders

	```sh
	cd ..
	mkdir select_seqs
	mkdir assembly
	```
2. Copy cleaned fastq files from KRLA/clean\_data to KRLA/select\_seqs

	```sh
	cd select_seqs
	nohup cp ../clean_data/*.fastq.gz . &> /dev/null &
	```
3. Check that all files have been moved correctly

	```sh
	ls *.fastq.gz -1 | wc -l
	```
	* **Number of KRLA individuals:** 497

### Generate unique sequence sets

There is no reason to use 100 identical sequences for the denovo clustering task, as it will only increase memory and energy usage and runtime.

1. Make a list of individual IDs from the .fastq.gz files

	```sh
	ls *.fastq.gz | sed -e 's/.fastq.gz//g' > namelist
	```

	* `-e` identifies the expression for sed
	* `'s/.fastq.gz//g'` substitute occurrances of '.fast.gz` with nothing, for all occurances

2. For each individual, create a file with only the unique reads from that individual (and a count of their occurances). This will speed up future steps.

	2a. Define variables to use with awk and perl:
	
	```sh
	AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
	AWK2='!/>/'
	PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
	```
	
	* The first two variables will be used in an awk command. They define the text processing commands that will be applied:
		* AWK1 defines a command that will look at every chunk of 4 lines (1 read of a fastq file). It will extract/print the first 2 lines of that chunk, and replace the beginning '@' with '>'
		* AWK2 defines a command that will ignore any lines that don't begin with '>'
	* The last variable defines a perl script that will be run later
		* The first portion loops through all the lines and stores/hashes them in `%z`, while incrementing a count for each value (effectively counting the number of times the sequence has occurred within the file)
		* The second portion iterates through the hash, and prints the count and then the sequence (tab separated)

	2b. Run the awk and perl commands:
	
	```sh
	nohup cat namelist | parallel --no-notice -j 8 "zcat {}.fastq | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs" &> /dev/null &
	```
	
	* This command gives the individual ID lines to the command (inside the "") run by parallel:
		* `parallel  --no-notice -j 8` runs the following command in parallel. `-j` is used to specify the number of jobs to run in parallel. `--no-notice` disables command line notices generated by `parallel`.
		* `zcat` is just like `cat` but for compressed files. It combines the individual ID with '.fastq.gz' to identify the file that contains that individual's reads
		* The following commands apply the awk and perl commands to the file. `mawk` is a version of `awk`, and `perl` runs the perl script.
		* The combination of commands takes the sequence lines from the .fastq.gz file, counts how many times they occur, and produced a count-sequence file that is saved in a file with the individual ID (stored in {}) with the .uniq.seqs ending
	* The command uses `nohup` and `&> /dev/null &` as described above

3. Check that you have a uniq.seq file for each fastq.gz file:

	```sh
	ls *.uniq.seqs | sed -e 's/.fastq.gz//g' > nameList
	```
	
### Assemble from sequence sets
	
4. Select a value of i (number of individuals a sequence occurs in) and k (number of times a sequence appears in an individual) to filter your .uniq.seqs files. This will speed up the assembly step.

	* Values are usually between 2 and 10
	* See note below on how to iteratively repeat steps 4-6 with different parameter values


5. Generate the collection of sequences that meet your i and k criteria by running the `selectContigs.sh` script. For example when k=4 and i=2:

	```sh
	cp /working/romero/scripts/selectContigs.sh .
	nohup bash selectContigs.sh 4 2 > k4.i2.seqs &> /dev/null &
	```
	
	* The file k4.12.seqs will contain only sequences that meet these criteria
	* selectContigs.sh contents:

		```sh
		#!/bin/bash
		
		parallel --no-notice -j 16 mawk -v x=$1 \''$1 >= x'\' ::: *.uniq.seqs \
    	| cut -f2 \
    	| perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' \
    	| mawk -v x=$2 '$1 >= x' \
    	| cut -f2 \
    	| mawk '{c= c + 1; print ">Contig_" c "\n" $1}' \
    	| sed -e 's/NNNNNNNNNN/\t/g' \
    	| cut -f1
		```

6. Use [CD-HIT](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#user-content-CDHITEST) to create a denovo assembly from these reads, at a chosen clustering similarity threshold (c).
	
	```sh
	module load cd-hit/4.6
	nohup cd-hit-est -i <inputFile> -o <outputFile> -M 0 -T 0 -c 0.92 &>/dev/null &
	```
	
	* `<inputFile>` is your file from step 5 (k4.i2.seqs)
	* `<outputFile>` is the filename you want for your output. CD-HIT will generate 2 output files:
		* 1) outputFile (no extension)
		* 2) outputFile.clstr
	* `-M` is max memory allowance, 0 sets it to unlimited, default is 800M
	* `-T` is max number of threads / CPUs, 0 sets it to unlimited (32), it may be better to set it to ~16 to not take up the entire server
	* `-c` is the clustering similarity (how similar sequences need to be to be joined together)
		* 0.95 is usually a good clustering similarity to use, but you should take into account how diverged the populations or species you are studying are
		* Higher values of c create more contigs and runs faster
		* Lower values of c create fewer contigs and runs slower

**Note 1:** If you are interested in comparing the results of differnt combinations of i, k, and c parameters effect the resultant number of contigs, you can use the genContigSets.sh script to create kn.in.seqs files for many combinations (each combination of i and k across 2,4,6,8, and 10), and then run each of those files through CD-HIT and compare results. However, it is difficult to interpret the number of contigs (if it is over- or under-assembled).

* If you do compare results from multiple combinations of parameters, it would make sense to set your <outputFile> name to include the i, k, and c parameters (like rf4.2.92 for k = 4, i = 2, and c = 0.92, this will also allow easy detection of assembly files to compare number of contigs)
* **Will add information later on a script that parallelizes multiple cd-hit assemblies for comparison...**
* To summarize the information from different assemblies:
	
	```sh
	grep "^>" rf*[0-9] -c | awk -F"[:.]" '{print $2"\t"$3"\t"$4"\t"$5}' > assemblyComparison
	less assemblyComparison
	```
	
	* **Will add some code and/or images of plots for these comparisons later**

**Note 2:** Another option for comparing assembly parameters is.....

* what the 'refOpt.sh' attempts to do (i.e. what Trevor does on pronghorn). 
* testing parameter effects across a subset of individuals vs all individuals
* Processing time (day-ish) and disk space is honestly not that big relative to this whole pipeline and pretty feasible on ponderosa (with potential to improve efficiency still)
* Might be worth doing in this case just to get a more in-depth understanding of things even if we decide it's not worth it on future projects...



# 3. Mapping reads from all individuals to reference using `bwa`

### 3A. DIRECTORY & FILE PREP

## III. Using bcftools to build cigar formatted mpileup

## IV. Estimating genotype likelihoods with SAMtools

## V. Filtering variants with VCFtools

## VI. Entropy
---
---
BELOW IS EXAMPLE OF TWO APPROACHES OF ALIGNING TO REFERENCE GENOME:
## Alignment to *T. cristinae* genome and variant calling.
New versions of software installed on ponderosa, with modules:
- bwa 0.7.17-r1188 (https://github.com/lh3/bwa/releases)
- bcftools 1.9 (under https://sourceforge.net/projects/samtools/files/samtools/1.9/)
- samtools 1.10 (under https://sourceforge.net/projects/samtools/files/samtools/1.10/)

## 1) Working with T. cristinae reference genome
located at:
    ponderosa:/working/parchman/tpodura/raw_ind_fastqs/
    
## 2) reference based assembly with `bwa 0.7.17-r1188`

`NOTE`: Moving forward with 598 individuals

Make index for `bwa`

    $ module load bwa/0.7.17-r1188
    $ bwa index -p cristinae -a bwtsw re_mod_map_timema_06Jun2016_RvNkF702.fasta &


`bwa` wrapper, runbwa_memTLP.pl, modified to run the `mem` algorithim (rather than aln and samse), and used bwa 0.7.17-r1188. Parameter settings are described within the wrapper, more info with `bwa mem`

    $ module load bwa/0.7.17-r1188
    $ perl runbwa_memTLP.pl  *fastq &

## 3) Sorting, indexing, and converting `sam` files to `bam`

Number of threads set in script, based on current server usage.

    $ module load samtools/1.10
    $ perl sam2bamV1.10.pl *.sam


Cleaning up the directory
    
    $ rm *.sam
    $ rm *.sai

## 4. Making pileup (bcf) and variant calling with `bcftools 1.9`
tpod_bams is a text file with all of the **sorted** bam files listed, one per line

Options used:

-C --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]

-d --max-depth INT     max per-file depth; avoids excessive memory usage [250] 

-f --fasta-ref FILE    faidx indexed reference sequence file

-q --min-MQ INT        skip alignments with mapQ smaller than INT [0]

-Q --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]

-I --skip-indels       do not perform indel calling

-b --bam-list FILE     list of input BAM filenames, one per line

-O --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;
                          'z' compressed VCF; 'v' uncompressed VCF [v]
                          
-o --output FILE       write output to FILE [standard output]


    $ module load bcftools/1.9
    $ bcftools mpileup -C 50 -d 250 -f re_mod_map_timema_06Jun2016_RvNkF702.fasta -q 30 -Q 20 -I -b tpod_bams -O b -o tpod.bcf

## 5. Generation vcf file from bcf.
Options used:
-v --variants-only             output variant sites only

-c --consensus-caller          the original calling method (conflicts with -m)

-f --format-fields <list>      output format fields: GQ,GP (lowercase allowed) []
    
-p --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5]
                                                         
-P --prior <float>         mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
    
-O --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

    $ module load bcftools/1.9
    $ bcftools call -v -c -f GQ -p 0.01 -P 0.001 -O v -o tpod.vcf tpod.bcf

Checking number of SNPs:
    
    $ module load vcftools/0.1.14
    $ vcftools --vcf tpod.vcf 

After filtering, kept 598 out of 598 Individuals
After filtering, kept 127722 out of a possible 127722 Sites

## 5. Filtering

Just doing some rough preliminary stuff here, need to consider how ZG recommends filtering based on what was done above, and how he has been doing things with Timema for mapping.

make id file for reheadering

    $ ls *fastq > fastqs.txt
    $ sed -s "s/.fastq//" fastqs.txt > Tpod_ids_col.txt

reheader vcf

    $ module load bcftools/1.9
    $ module load vcftools/0.1.14
    $ bcftools reheader -s Tpod_ids_col.txt tpod.vcf -o rehead_tpod.vcf

initial round of filtering (just getting a feeling of how filtering parameters might shape the dataset)

    $ vcftools --vcf rehead_tpod.vcf --out variants_maf5_miss5 --remove-filtered-all --maf 0.03 --max-missing 0.5 --recode --thin 100

After filtering, kept 598 out of 598 Individuals, 18640 out of a possible 127722 Sites

    $ vcftools --vcf rehead_tpod.vcf --out variants_maf3_miss5 --remove-filtered-all --maf 0.03 --max-missing 0.5 --recode --thin 100
	        

After filtering, kept 598 out of 598 Individuals, kept 19384 out of a possible 127722 Sites


## Due to missing 702.1 in *T. cristinae* genome, now trying alignment to *T. podura* consensus genome and variant calling.

## 5) Working with *T. podura* reference genome
located at:
    ponderosa:/working/parchman/tpodura/raw_ind_fastqs/Tpodura_consensus.fa
    
## 6) reference based assembly with `bwa 0.7.17-r1188`

`NOTE`: Moving forward with 598 individuals

Make index for `bwa`

    $ module load bwa/0.7.17-r1188
    $ bwa index -p Tpodura_consensus -a bwtsw Tpodura_consensus.fa &


`bwa` wrapper, runbwa_memTLP.pl, modified to run the `mem` algorithim (rather than aln and samse), and used bwa 0.7.17-r1188. Parameter settings are described within the wrapper, more info with `bwa mem`

    $ module load bwa/0.7.17-r1188
    $ perl runbwa_memTLP.pl  *fastq &

## 7) Sorting, indexing, and converting `sam` files to `bam`

Number of threads set in script, based on current server usage.

    $ module load samtools/1.10
    $ perl sam2bamV1.10.pl *.sam


Cleaning up the directory
    
    $ rm *.sam
    $ rm *.sai

## 8. Making pileup (bcf) and variant calling with `bcftools 1.9`
tpod_bams is a text file with all of the **sorted** bam files listed, one per line

Options used:

-C --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]

-d --max-depth INT     max per-file depth; avoids excessive memory usage [250] 

-f --fasta-ref FILE    faidx indexed reference sequence file

-q --min-MQ INT        skip alignments with mapQ smaller than INT [0]

-Q --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]

-I --skip-indels       do not perform indel calling

-b --bam-list FILE     list of input BAM filenames, one per line

-O --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;
                          'z' compressed VCF; 'v' uncompressed VCF [v]
                          
-o --output FILE       write output to FILE [standard output]


    $ module load bcftools/1.9
    $ bcftools mpileup -C 50 -d 250 -f Tpodura_consensus.fa -q 30 -Q 20 -I -b tpod_bams -O b -o tpod.bcf

# Done to here.

## 9. Generation vcf file from bcf.
Options used:
-v --variants-only             output variant sites only

-c --consensus-caller          the original calling method (conflicts with -m)

-f --format-fields <list>      output format fields: GQ,GP (lowercase allowed) []
    
-p --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5]
                                                         
-P --prior <float>         mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
    
-O --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

    $ module load bcftools/1.9
    $ bcftools call -v -c -f GQ -p 0.01 -P 0.001 -O v -o tpod.vcf tpod.bcf

Checking number of SNPs:
    
    $ module load vcftools/0.1.14
    $ vcftools --vcf tpod.vcf 

After filtering, kept 598 out of 598 Individuals
After filtering, kept 127722 out of a possible 127722 Sites

## 10. Filtering

Just doing some rough preliminary stuff here, need to consider how ZG recommends filtering based on what was done above, and how he has been doing things with Timema for mapping.

make id file for reheadering

    $ ls *fastq > fastqs.txt
    $ sed -s "s/.fastq//" fastqs.txt > Tpod_ids_col.txt

reheader vcf

    $ module load bcftools/1.9
    $ module load vcftools/0.1.14
    $ bcftools reheader -s Tpod_ids_col.txt tpod.vcf -o rehead_tpod.vcf

initial round of filtering (just getting a feeling of how filtering parameters might shape the dataset)

    $ vcftools --vcf rehead_tpod.vcf --out variants_maf5_miss5 --remove-filtered-all --maf 0.03 --max-missing 0.5 --recode --thin 100

After filtering, kept 598 out of 598 Individuals, 18640 out of a possible 127722 Sites

    $ vcftools --vcf rehead_tpod.vcf --out variants_maf3_miss5 --remove-filtered-all --maf 0.03 --max-missing 0.5 --recode --thin 100
	        

After filtering, kept 598 out of 598 Individuals, kept 19384 out of a possible 127722 Sites

### Appendix 1: Reference of useful commands (as far as what Parchman lab people like)
Just an unorganized list for now, will clean up later...

+ `parallel` - easiest way to run jobs in parallel across CPUs. See [GNU Parallel tutorial](https://www.gnu.org/software/parallel/parallel_tutorial.html)
    + `-j` - max jobs to run in parallel
    + `--no-notice` - eliminates excessive notices printing to shell when running jobs
    + `:::` - followed by list of input files to process
+ `nohup <command> &> /dev/null &` - a way to run a longer background process that won't be interupted
    + `nohup` - keeps process running even if shell or terminal is exited (i.e. long jobs don't get terminated)
    + `&` - process is put in background (have access to shell while process is running)
    +  `/dev/null` - essentially a black hole to direct st. output from nohup into assuming you've already captured the output of interest
    + can do similar things with `screen` but `nohup` is simpler and enough for most of the use cases here
+ `time <command>` - prints the time a process takes after it completes
    + will generate 3 different times, but "real" is what you're usually interested in
    + useful for testing pararmeters of parallelization and getting idea of how long different tasks in pipeline take
+ `du -sch <files/dir/etc.> | tail -n 1` - way to see how much disk space a set of files is using, useful if a lot of temporary/intermediate files are being generated
+ `htop` - monitor status of jobs and CPU usage (just google for details)
