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