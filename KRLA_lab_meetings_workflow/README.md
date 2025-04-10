# Genotyping By Sequencing (GBS) Example Study

This document presents our workflow and rationale for genotype inference from high throughput sequencing of reduced representation libraries (GBS, RADseq, ddRADseq, etc). Several canned software packages or computational workflows exist for handling this type of data. These methods rely on set thresholds for sequencing coverage depth per locus to call hard genotypes. The biggest cost of using these methods is throwing away much, if not most, of the data. Examples of such software/workflows include:

* [Stacks](http://catchenlab.life.illinois.edu/stacks/)
* [Ddocent](http://www.ddocent.com//)
* [ipyrad](https://ipyrad.readthedocs.io/en/master/)

See [Nielsen et al. 2011](./papers/Nielsen_etal_2011.pdf) and [Buerkle and Gompert 2013](./papers/Buerkle_Gompert_2013.pdf) for articulate thoughts about this.

## Running to do list

* insert plate map and ID file path in DNA extractions
* add upstream path of mastermix cocktail file
* fill in number of reads before and after cleaning
* Add explanation of selectContigs.sh parts
* Add estimated run times for each step

## Table of Contents

 1. Species Information
 2. Directory Structure
 3. DNA Extraction
 4. Library Preparation
 5. Sequencing
 6. Data Cleaning
    1. Cleaning Contaminants
    2. Barcode Parsing
    3. Splitting fastqs
 7. Denovo Assembly
    1. Prepare Directories and Files
    2. Generate Unique Sequence Sets
    3. Assemble from Sequence Sets
 8. Mapping Reads with BWA
    1. Prepare Directories and Files
    2. Map, Sort, and Index
    3. Explanation of Script
 9. Build Mpileup with bcftools
    1. Prepare Directories
    2. Pileup, Call, and Filter
    3. Understanding bcftools Parameters
10. Convert bcf to vcf
    1. Identify SNPs
    2. Understanfing vcftools Parameters
11. Filtering
    1. Filtering on Individuals
    2. Filtering on Loci
12. Uncategorized
    1. Reference-based assembly (*T. podura*)
    2. Reheader vcf files?
13. Appendices
    1. Appendix 1: Useful Commands

## Focal Species: *Krascheninnikovia lanata*

<img src="images/KRLAplant.jpg" alt="A Krascheninnikovia lanata plant" width="400">
<img src="images/KRLAmap.png" alt="Map of sampling locations" width="400">

*Krascheninnikovia lanata* (winterfat) is a perennial shrub with a broad north/south distributional range spanning western Canada, U.S. and Mexico. The species is exclusive to North America and its current range is likely the result of southward expansion following two distinct migration events from eastern Mongolian lineages ~ 1.8 - 0.5 Mya.

The species is a halophyte (salt-tolerant) and is one of the only species outside of the *Atriplex* complex to co-dominate the salt desert shrublands of the Great Basin. It is a highly nutritious source of forage which is notable given that it is prone to replacement by the toxic exotic *Halogeton glomeratus* within disturbed habitats. The common name 'winterfat' is indicative of the persistence of green leaves throughout the winter season and late fall phenology.

Population sampling was a combined effort throughout 2021 - 2022 with Cathy Silliman doing collections for most of the populations in the west, central, and north Great Basin and Seth Romero gathering collections from the eastern Great Basin and Mojave. Anecdotally, many of the individuals in the north, central and eastern Great Basin were smaller in stature but part of broad near-monocultures that created consistent cover across broad areas. By contrast, many of the populations in the Mojave and southwest Great Basin were large individuals that occured in small islands or as sub-dominants with only a handful of individuals living in close proximity. The populations **DT** and **CL**, in particular, had nearly every individual sampled that could found at those locations.

## Directory Structure

* Create a project folder in your working directory called KRLA with `mkdir KRLA`
* The following directory structure should be made throughout the tutorial

```mermaid
flowchart TD;
    A(personal directory <br> /working/romero/) --> B(species folder <br> /romero/KRLA/)
    B --> C(assembly)
    B --> N (select_seqs)
    B --> D(bwa)
    B --> E(fastq)
    B --> F(scripts)
    E --> G(fastq files <br> e.g. *.fastq.gz)
    C --> H(de novo assembly <br> e.g. rf.*)
    C --> I(indexed assembly files <br> e.g. *.amb, *.ann, *.bwt, *.pac, *.sa)
    C --> J(alt_assemblies <br> OPTIONAL)
    J --> K(seq subset files <br> e.g. k4.i2.seqs)
    J --> L(assembly files <br> e.g. rf.4.2.92)
    D --> M(mapped/sorted reads <br> + index files <br> e.g. *.bam and *.bam.bai)
```

## DNA extraction

* Performed by AG Biotech in 12/2022
* DNA information with plate maps and IDs can be found at [insert doc path here]

## Library preparation

* Performed in Parchman Lab in 12/2022
* [add upstream path]KRLA\_RFseq\_mastermixcockatils.xlsx contains information about reagents used in library prep
* R/L and PCR for plates 1-6

## Sequencing

* 1 lane of S2 chemistry NovaSeq data at UTGSAF in 03/2023
* sequence results in /archive/parchman\_lab/rawdata\_to\_backup/KRLA/KRLA\_S1\_L001\_R1\_001.fastq.gz

## Data cleaning

### Cleaning contaminants

* **Goal:** Remove reads that are not from the target organims
* Performed 06/09/23

1. Create a fastq folder inside KRLA folder:

   ```sh
   mkdir fastq
   ```

2. Copy the sequencing results to your fastq folder:

   ```sh
   cd fastq
   cp /archive/parchman_lab/rawdata_to_backup/KRLA/KRLA_S1_L001_R1_001.fastq.gz .
   ```

3. Decompress the file with:

   ```sh
   gunzip KRLA_S1_L001_R1_001.fastq.gz
   ```

4. Count the reads **before** cleaning with:

   ```sh
   nohup grep -c "^@" KRLA_S1_L001_R1_001.fastq > KRLA_number_of_rawreads.txt &
   ```

   * `nohup` (no hang up) allows the command to keep running, even if you close the terminal, lose connection with the server, or log out.
   * `&` allows the command to run in the background.
   * `grep -c` will count the number of occurances of the pattern ("^@") in the file (KRLA_S1_L001_R1_001.fastq).
   * the result file (KRLA_number_of_rawreads.txt) will contain the initial number of reads.
   * **Reads before cleaning:** n

5. Run the cleaning script (note: the script will need to be edited to change the paths to the appropriate files):

   ```sh
   nano ../scripts/cleaning_bash_KRLA.sh # this is where you edit the path names
   conda deactivate # to make sure you are using the loaded modules
   module load fqutils/0.4.1
   module load bowtie2/2.2.5
   bash cleaning_bash_KRLA.sh &
   ```

   * Run this in the fastq directory
   * `module load` is needed to allow access to the modules used in the script
   * fqutils and bowtie2 are used in the cleaning script
   * cleaning\_bash\_KRLA.sh uses the tapioca pipeline's contamination analysis script to remove artifacts from:
      * Illumina oligos
      * PhiX
      * *E. coli*

6. Check that KRLA.clean.fastq has been created (from cleaning script), then remove the raw data file (from KRLA/fastq - a backup is stored at /archive/parchman\_lab/rawdata\_to\_backup/KRLA) :

   ```sh
   ls # look for KRLA.clean.fastq
   rm -rf KRLA_S1_L001_R1_001.fastq &
   ```

7. Count the reads **after** cleaning with:

   ```sh
   nohup cd-hit-est -i k4.i2.seqs -o rf.4.2.92 -M 0 -T 0 -c 0.92 &>/dev/null &
   ```

   ```sh
   nohup grep -c "^@" KRLA.clean.fastq > KRLA_clean_reads.txt &
   ```

   * reads after cleaning are typically ~75% of reads before cleaning
   * **Reads after cleaning:** n

### Barcode parsing

* **Goal:** Remove barcode and adapter sequences from reads, place this information (what individual the reads came from) in the read ID line

1. Run parsing script:

   ```sh
   cp /working/parchman/KRLA/parse_barcodes769.pl ../scripts
   cp /working/parchman/KRLA/KRLA_barcode_key.csv ../scripts
   nohup perl ../scripts/parse_barcodes769.pl ../scripts/KRLA_barcode_key.csv KRLA.clean.fastq A00 &>/dev/null &
   ```

   * Run this in the fastq directory
   * KRLA\_barcode\_key.csv provides the barcodes for each individual
   * KRLA.clean.fastq are the clean reads from contaminant cleaning
   * A00 is the sequencer ID (first 3 characters after @ in the fastq identifier)
   * &>/dev/null prevents display of stdout

2. View the parsing results (created by the parsing script):

   ```sh
   less parsereport_KRLA.clean.fastq
   ```

   * Results compare good and bad mids (molecular IDs), bad mids are usually <5% of total mids
   * The other removed results are typically insignificant
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
   cp /working/parchman/KRLA/splitfastqs/splitFastq_universal_regex.pl ../scripts
   nohup perl ../scripts/splitFastq_universal_regex.pl KRLA_ids_noheader.txt parsed_KRLA.clean.fastq &>/dev/null &
   ```

   * Run this in the fastq directory
   * The perl script creates a separate fastq file for each individual's reads

3. Compress all of the resultant fastq files:

   ```sh
   nohup gzip *fastq &>/dev/null &
   ```

## Denovo assembly

Denovo assemblies uses your reads to generate a consensus artificial reference to map your reads against. This is the only option in some systems where there is no reference genome available. Sometimes it is still prefered, depending on the genome divergence between your species or population and that of the reference genome. If you are using a reference genome, download the reference and skip this step.

There are many tools available for denovo assembly. We will be using [CD-HIT](https://www.bioinformatics.org/cd-hit/) in a workflow from [dDocent](http://www.ddocent.com/) (Puritz et al. 2014). See [LaCava et al. 2020](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13108) for a comparison of tools.

This workflow begins with the gziped .fastq files in KRLA/fastq.

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
   nohup cp ../fastq/*.fastq.gz . &> /dev/null &
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

   * Run this in select_seqs
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

   * Run this in select_seqs
   * This command gives the individual ID lines to the command (inside the "") run by parallel:
      * `parallel  --no-notice -j 8` runs the following command in parallel. `-j` is used to specify the number of jobs to run in parallel. `--no-notice` disables command line notices generated by `parallel`.
      * `zcat` is just like `cat` but for compressed files. It combines the individual ID with '.fastq.gz' to identify the file that contains that individual's reads
      * The following commands apply the awk and perl commands to the file. `mawk` is a version of `awk`, and `perl` runs the perl script.
      * The combination of commands takes the sequence lines from the .fastq.gz file, counts how many times they occur, and produced a count-sequence file that is saved in a file with the individual ID (stored in {}) with the .uniq.seqs ending
   * The command uses `nohup` and `&> /dev/null &` as described above

3. Check that you have a uniq.seq file for each fastq.gz file:

   ```sh
   ls *.uniq.seqs -1 | wc -l
   ```

### Assemble from sequence sets

1. Select a value of i (number of individuals a sequence occurs in) and k (number of times a sequence appears in an individual) to filter your .uniq.seqs files. This will speed up the assembly step.

   * Values are usually between 2 and 10
   * See note below on how to iteratively repeat steps 1-3 with different parameter values

2. Generate the collection of sequences that meet your i and k criteria by running the `selectContigs.sh` script. For example when k=4 and i=2:

   ```sh
   cp /working/romero/scripts/selectContigs.sh ../scripts
   nohup bash ../scripts/selectContigs.sh 4 2 > k4.i2.seqs &
   ```

   * Run this in select_seqs
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

       * [insert description of what this is doing]

3. Use [CD-HIT](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#user-content-CDHITEST) to create a denovo assembly from these reads, at a chosen clustering similarity threshold (c).

   ```sh
   cd ../assembly
   module load cd-hit/4.6
   ```

   ```sh
   nohup cd-hit-est -i <inputFile> -o <outputFile> -M 0 -T 0 -c 0.94 &>/dev/null &
   ```

   * `<inputFile>` is your file from step 5 (../select_seqs/k4.i2.seqs)
   * `<outputFile>` is the filename you want for your output. CD-HIT will generate 2 output files:
      * 1) outputFile (no extension)
      * 2) outputFile.clstr
   * `-M` is max memory allowance, 0 sets it to unlimited, default is 800M
   * `-T` is max number of threads / CPUs, 0 sets it to unlimited (32), it may be better to set it to ~16 to not take up the entire server
   * `-c` is the clustering similarity (how similar sequences need to be to be joined together)
      * 0.95 is usually a good clustering similarity to use, but you should take into account how diverged the populations or species you are studying are
      * Higher values of c create more contigs and runs faster
      * Lower values of c create fewer contigs and runs slower

4. Use [bwa index](https://bio-bwa.sourceforge.net/bwa.shtml) to index our assembly into fasta format for mapping individual reads

   ```sh
   module load bwa/0.7.17-r1188
   bwa index -p K_lanata rf.4.6.94
   ```

   * `-p` lets us set the prefix name that will be attached to all output assembly files. Use something descriptive of the species you are working on.
   * `-a` (optional) allows you to choose the BWT construction algorithm
      * `is` (default) usually faster but limited to databases smaller than 2GB
      * `bwtsw` for use on large databases (e.g. human genome)
   * rf.4.6.94 is the denovo assembly from CD-HIT

   The result of this step should produce 5 new files named with your chosen prefix and the following extensions: `*.amb`, `*.ann`, `*.bwt`, `.pac`, `*.sa`

* **Note:** Above we describe the process for selecting one value of i, one value of k, and one value of c. If instead you would like to compare multiple assembly options, you can:

   1. Use the genContigSets.sh script to create kn.in.seqs files for many combinations, and then run each of those files through CD-HIT and compare results. However, it is difficult to interpret the number of contigs (if it is over- or under-assembled).
   2. Use refOpt.sh (what Trevor does on pronghorn), which:

      * tests parameter effects across a subset of individuals vs all individuals
      * takes about a day of processing time
      * make be feasible on ponderosa (disk space is not that big relative to this whole pipeline)
      * could give you a more in-depth understanding of things even if we decide it's not worth it on future projects

* If you do compare results from multiple combinations of parameters, it would make sense to set your <outputFile> name to include the i, k, and c parameters (like rf4.2.92 for k = 4, i = 2, and c = 0.92, this will also allow easy detection of assembly files to compare number of contigs)
* To summarize the information from different assemblies:

   ```sh
   grep "^>" rf*[0-9] -c | awk -F"[:.]" '{print $2"\t"$3"\t"$4"\t"$5}' > assemblyComparison
   less assemblyComparison
   ```

## Mapping reads with `bwa`

### Prepare Directories and Files

1. Make a new directory from the species base (e.g. `.../KRLA/`) called `bwa` and go into it

   ```sh
   cd ..
   mkdir bwa
   cd bwa
   ```

2. Move **indexed** assembly files into bwa directory. In this case these all have the prefix **"K_lanata"** based on the prior step.

   ```sh
   cp ../assembly/K_lanata* .
   ```

### Map, sort and index

1. Load both `bwa` and `samtools` which are required for running the mapping script

   ```sh
   module load bwa/0.7.17-r1188
   module load samtools/1.10
   ```

2. Run the script `bwa_mem2sorted_bam.sh` using the following nohup settings

   ```sh
   nohup bash ../scripts/bwa_mem2sorted_bam.sh 2> /dev/null &
   ```

   * Running the script in this way prevents the process from being interrupted (i.e. you can disconnect from the server while this runs) while also capturing progress print statements in `nohup.out`. You can re-login to the server and check the progress of mapping by going into the `bwa` directory and entering the following:

      ```sh
      tail -n 1 nohup.out
      ```

   * This step took **~6 hours** using **24 nodes** on ponderosa for **497 individuals** in the KRLA dataset.

### Explanation of `bwa_mem2sorted_bam.sh`

The contents of the previous script is the following:

```sh
#!/bin/bash

ctr=1
fTotal=$(ls ../fastq/*.gz -1 | wc -l)

for file in ../fastq/*.gz
   do
   if test -f "$file"
   then
      fPrefix=$(echo "$file" | sed 's|.*/||' | cut -f1 -d ".")
      echo mapping and sorting individual "$ctr" of "$fTotal"
      bwa mem -t 24 K_lanata "$file" | \
      samtools view -u -b - | \
      samtools sort -l0 -@24 -o "$fPrefix.sorted.bam"
      ((ctr++))
   fi
done
for sBam in *.sorted.bam
   do
   if test -f "$sBam"
   then
      samtools index -@24 "$sBam"
   fi
done
```

We should explain the steps that are happening here particularly any settings used with

* `bwa mem` - maps sequences to the reference, creating a .sam file
  * `-t` - number of threads used
* `samtools view` - converts .sam format to .bam format to save space
  * `-u` - output uncompressed data
  * `-b` - output in .bam format
* `samtools sort` - sorts .bam files by position on the reference
  * `-l` - output compression level (l0) = uncompressed
  * `-@` - number of threads used (@24 = 24 threads)
  * `-o` - output file name
* `samtools index` - indexes the .bam file for faster search, creating .bai files
  * `-@` - number of threads used (@24 = 24 threads)

## Build pileup and variant call with BCFTools

### Prepare Directoties

1. Make vcf directory

   ```sh
   cd ..
   mkdir vcf
   ```

2. Make list of bam files from

```sh
find /bwa/ -type f -name *.sorted.bam > bam_list.txt
```

### Pileup, call, and filter

1. Option 1 using `bcftools 1.3` (from Lainie)

   * The following takes several (2-4) hours. NOTE: Run the following lines as one large chunk of code. Be sure to change reference name and output file.

      ```sh
      samtools mpileup -P ILLUMINA --BCF --max-depth 130 --adjust-MQ 50 --min-BQ 20 --min-MQ 20 --skip-indels --output-tags DP,AD --fasta-ref K_lanata.fa *sorted.bam | \
      bcftools call -m --variants-only --format-fields GQ --skip-variants indels | \
      bcftools filter --set-GTs . -i 'QUAL > 19 && FMT/GQ >9' | \
      bcftools view -m 2 -M 2 -v snps --apply-filter "PASS" --output-type v  --output-file variants_rawfiltered_1FEB2024.vcf &
      ```

<<<<<<< Updated upstream
2. Option 2 using `bcftools 1.9`:

   ```sh
   module load bcftools/1.9
   bcftools mpileup -C 50 -d 250 -f K_lanata.fa -q 30 -Q 20 -I -b bam_list.txt -O b -o K_lanata.bcf
   ```
=======
Alternatively, we could use `bcftools 1.9` and run the following:
```sh
module load bcftools/1.9
```

```sh
bcftools mpileup -a DP,AD,INFO/AD -C 50 -d 250 -f K_lanata.fa -q 30 -Q 20 -I -b bam_list.txt -o KRLA.bcf
```

```sh
bcftools call -v -m -f GQ KRLA.bcf -O z -o KRLA.vcf.gz
```

```sh
module load vcftools/0.1.14
```

```sh
vcftools \
--remove-indels \
--min-alleles 2 \
--max-alleles 2 \
--thin 100 \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--gzvcf KRLA.vcf.gz \
--out 
```

```sh
awk '$5 > 0.5 {print $1}' KRLA.imiss | tail -n +2 > indmiss50.txt
```

module load bcftools/1.9
    $ bcftools mpileup -C 50 -d 250 -f K_lanata.fa -q 30 -Q 20 -I -b bam_list.txt -O b -o K_lanata.bcf
>>>>>>> Stashed changes

### Understanding bcftools parameters

* -C --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]
* -d --max-depth INT     max per-file depth; avoids excessive memory usage [250]
* -f --fasta-ref FILE    faidx indexed reference sequence file
* -q --min-MQ INT        skip alignments with mapQ smaller than INT [0]
* -Q --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
* -I --skip-indels       do not perform indel calling
* -b --bam-list FILE     list of input BAM filenames, one per line
* -O --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
* -o --output FILE       write output to FILE [standard output]

## Convert bcf file to vcf file

### Identify SNPs

1. Load Modules and Convert

   ```sh
   module load bcftools/1.9
   bcftools call -v -c -f GQ -p 0.01 -P 0.001 -O v -o tpod.vcf tpod.bcf
   ```

2. Check number of SNPs

   ```sh
   module load vcftools/0.1.14
   vcftools --vcf tpod.vcf 
   ```

   * After filtering, kept 598 out of 598 Individuals
   * After filtering, kept 127722 out of a possible 127722 Sites

### Understanding vcftools Parameters

* -v --variants-only             output variant sites only
* -c --consensus-caller          the original calling method (conflicts with -m)
* -f --format-fields <list>      output format fields: GQ,GP (lowercase allowed) []
* -p --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5]
* -P --prior <float>         mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
* -O --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

## Filtering

### Filtering on Individuals

* **Coverage:** Also can be thought of as depth. See 3mapping. Calculated on bam files. Average read count per locus per individual.
* **Missing:** Proportion of missing data allowed across all loci for individual. Common and high in GBS/RADseq data. Kinda an issue all around. Many methods, including PCA (all ordination methods), require a complete matrix with no missing data. Additionally, PCA will cluster by missing data with individuals with higher missing data clustering closer to the center and get this "fan" effect. Can be the same for coverage too. This (among other reasons) is why people use a variance-covariance matrix of genetic data to do ordinations. Other methods involve imputation. This can be fancy and use phased haplotype data OR simply, when you z-score, (g - mean(g))/sd(g), your genotype data across each locus, you make all missing data equal to 0 or Mean (i.e., the global allele frequency). There's more to this standardization, see [Patterson et al. 2006](https://dx.plos.org/10.1371/journal.pgen.0020190) for more info. See PCAsim_ex in examples directory for showing all these issues. This is another reason to use entropy. Entropy is a hierarchical bayesian model so it gets an updated genotype estimate for each missing value based on genotype likelihoods across loci, individuals, and the allele frequency of the cluster/deme that individual assigns to.

### Filtering on Loci

* **Biallelic:** Only keep biallelic SNPs. Multiallelic SNPs are rare at the time scale we work (Citation??) and also, mathematical nightmare and we have enough data so just ignore. Everyone does unless deep time phylogenetics.
* **thin:** Keeps one locus within a specified range. Not 100% how it decides with one to keep. I think it's on quality or depth. This is a necessary step as loci in close physical are prone to sequencing error and linkage disequalibrium (LD) confounds many different population genetic parameters. For de novo reference assemblies, we thin to 100 as contigs/reads are ~92 bp in length. This keeps one locus per contig to control for LD and sequencing error, which is really common in pop gen and necessary for many analyses.
* **max-missing** = max proportion of missing data per locus
* **MAF** = minor allele frequency. Proportion of individuals a alternate allele needs to be present in order for that locus to be kept as a SNP. (e.g. maf = 0.02 for 250 individuals means that an alternate allele needs to be present in at least 5 individuals to be kept) Many papers have shown this is a big issue in clustering and demography (Citation). We do this a second time near the end if we removed individuals during missing data filtering.
* **Mean Depth:** Average allelic depth or read depth per locus. Too low could be sequencing error, too high could be PCR replication artifact (Citation).
* **Qual:** Locus Quality. You can look up the math. Usually above 20-30 is good but given our coverage and number of individuals, we can usually go way higher.
* **Fis:** Inbreeding coefficient. This is a contentous topic. This has to do with paralogs or paralogous loci. This is where loci map to multiple regions of the genome. Issues in highly repeative genomes. Usually leads to an excess of heterozygotes. Filtering on negative Fis can help. See these two McKinney papers [1](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12763), and [2](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12613). Katie and others in the lab use his package called HDPlot to deal with this.

## Uncategorized

### Reference-based assembly (*T. podura*)

1. Refernce file at: /working/parchman/tpodura/raw_ind_fastqs/Tpodura_consensus.fa
2. Make index for `bwa`

   ```sh
   module load bwa/0.7.17-r1188
   bwa index -p Tpodura_consensus -a bwtsw Tpodura_consensus.fa &
   ```

    * `bwa` wrapper, runbwa_memTLP.pl, modified to run the `mem` algorithim (rather than aln and samse), and used bwa 0.7.17-r1188. Parameter settings are described within the wrapper, more info with `bwa mem`

      ```sh
      module load bwa/0.7.17-r1188
      perl runbwa_memTLP.pl  *fastq &
      ```

### Reheader vcf files?

1. Make id file for reheadering

   ```sh
   ls *fastq > fastqs.txt
   sed -s "s/.fastq//" fastqs.txt > Tpod_ids_col.txt
   ```

2. Reheader vcf

   ```sh
   module load bcftools/1.9
   module load vcftools/0.1.14
   bcftools reheader -s Tpod_ids_col.txt tpod.vcf -o rehead_tpod.vcf
   ```

3. Initial filtering

   ```sh
   vcftools --vcf rehead_tpod.vcf --out variants_maf5_miss5 --remove-filtered-all --maf 0.03 --max-missing 0.5 --recode --thin 100
   ```

## Appendices

### Appendix 1: Reference of useful commands (as far as what Parchman lab people like)

* `parallel` - easiest way to run jobs in parallel across CPUs. See [GNU Parallel tutorial](https://www.gnu.org/software/parallel/parallel_tutorial.html)
* `-j` - max jobs to run in parallel
* `--no-notice` - eliminates excessive notices printing to shell when running jobs
* `:::` - followed by list of input files to process
* `nohup <command> &> /dev/null &` - a way to run a longer background process that won't be interupted
  * `nohup` - keeps process running even if shell or terminal is exited (i.e. long jobs don't get terminated)
  * `&` - process is put in background (have access to shell while process is running)
  * `/dev/null` - essentially a black hole to direct st. output from nohup into assuming you've already captured the output of interest
  * can do similar things with `screen` but `nohup` is simpler and enough for most of the use cases here
* `time <command>` - prints the time a process takes after it completes
  * will generate 3 different times, but "real" is what you're usually interested in
  * useful for testing pararmeters of parallelization and getting idea of how long different tasks in pipeline take
* `du -sch <files/dir/etc.> | tail -n 1` - way to see how much disk space a set of files is using, useful if a lot of temporary/intermediate files are being generated
* `htop` - monitor status of jobs and CPU usage (just google for details)
