# Organization and Workflow for *Krascheninnikovia lanata* GBS 
Organizational notes and code for =rangewide sampling for landscape genomic analyses

# Range-wide landscape genomics: sample organization and GBS workflow 

## Sample organization
- Full information on DNAs for each individual sampled across natural distribution can be found in `XXXXXXXXXX`. This file also has the updated plate maps with specified IDs.

- **NOTE** DNA was extracted in December 2022 at AG Biotech. Plates in lab freezer need to be tranported to -80.

## Notes on library preparation

### 12/19-12/22: R/L and PCR for plates 1-6. Master mix in `KRLA_RFseq_mastermixcockatils.xlsx`.


## Data analysis: contaminant cleaning, barcode parsing, data storage, directory organization, and initial analyses.

We generated 1 lane of S2 chemistry NovaSeq data at UTGSAF in March of 2023. 

## This file contains code and notes for
1) cleaning contaminants using tapioca
2) parsing barcodes
3) splitting fastqs 
4) de novo assembly
5) reference based assembly
6) calling variants
7) filtering
8) entropy for genotype probabilities.

## 1. Cleaning contaminants

Being executed on ponderosa using tapioca pipeline. Commands in two bash scripts (cleaning_bash_CADE.sh and cleaning_bash_SEGI.sh), executed as below (6/9/23). This was for one S2 NovaSeq lanes generated in late December 2022.

Decompress fastq file:

    $ gunzip KRLA_S1_L001_R1_001.fastq.gz

Number of reads **before** cleaning:

    $ nohup grep -c "^@" KRLA_S1_L001_R1_001.fastq > KRLA_number_of_rawreads.txt &
    ## raw reads: 

To run cleaning_bash* tapioca wrapper, exit conda environment, load modules, and run bash scripts.

    $ module load fqutils/0.4.1
    $ module load bowtie2/2.2.5
    
    $ bash cleaning_bash_KRLA.sh &


After .clean.fastq has been produced, rm raw data:

    $ rm -rf KRLA_S1_L001_R1_001.fastq &
   

Raw data will stay stored in: /archive/parchman_lab/rawdata_to_backup/FRLA/

Number of reads **after** cleaning:

    $ nohup grep -c "^@" KRLA.clean.fastq > FRLA1_clean_reads.txt &
    ## reads after cleaning:


## 2. Barcode parsing:


Be sure to deactivate conda environment before running the below steps. Barcode keyfiles are `/working/parchman/KRLA/KRLA_barcode_key.csv` 

Parsing KRLA library:

    $ nohup perl parse_barcodes768.pl KRLA_barcode_key.csv KRLA.clean.fastq A00 &>/dev/null &




`NOTE`: the A00 object is the code that identifies the sequencer (first three characters after the @ in the fastq identifier).

    $ less parsereport_KRLA.clean.fastq
    Good mids count: 1571963061
    Bad mids count: 73508410
    Number of seqs with potential MSE adapter in seq: 321195
    Seqs that were too short after removing MSE and beyond: 428

## 3. splitting fastqs


For KRLA, doing this in `/working/parchman/KRLA/splitfastqs`

Make ids file

    $ cut -f 3 -d "," KRLA_barcode_key.csv | grep "_" > KRLA_ids_noheader.txt


Split fastqs by individual

    $ nohup perl splitFastq_universal_regex.pl KRLA_ids_noheader.txt parsed_KRLA.clean.fastq &>/dev/null &


'gzip' all .fastq files

gzipped the parsed*fastq files for now, but delete once patterns and qc are verified.

    $ nohup gzip *fastq &>/dev/null &
    

# Workflow for assembly and variant calling for GBS data
Here we will organize thoughts, details, justification, and code for all steps, from raw data to final filtered genotype matrices for population genetic analyes

## Assembly

- when to use reference based
- when to not use reference based and why
- how to run reference based assembly
- how to run denovo assembly to generate artificial reference

####################################################################
####################################################################
BELOW IS EXAMPLE OF TWO APPROACHES OF ALIGNING TO REFERENCE GENOME:
## Alignment to *T. cristinae* genome and variant calling.
New versions of software installed on ponderosa, with modules:
- bwa 0.7.17-r1188 (https://github.com/lh3/bwa/releases)
- bcftools 1.9 (under https://sourceforge.net/projects/samtools/files/samtools/1.9/)
- samtools 1.10 (under https://sourceforge.net/projects/samtools/files/samtools/1.10/)

## 5) Working with T. cristinae reference genome
located at:
    ponderosa:/working/parchman/tpodura/raw_ind_fastqs/
    
## 6) reference based assembly with `bwa 0.7.17-r1188`

`NOTE`: Moving forward with 598 individuals

Make index for `bwa`

    $ module load bwa/0.7.17-r1188
    $ bwa index -p cristinae -a bwtsw re_mod_map_timema_06Jun2016_RvNkF702.fasta &


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
    $ bcftools mpileup -C 50 -d 250 -f re_mod_map_timema_06Jun2016_RvNkF702.fasta -q 30 -Q 20 -I -b tpod_bams -O b -o tpod.bcf

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