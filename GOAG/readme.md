## Desert tortoise analysis notes

### Notes on contaminant cleaning and barcode parsing 10-30

`NOTE`: One library was sequenced on one Novaseq lane in late September of 2021
`NOTE`: Contaminant cleaning and barcode parsing in `/working/parchman/GOAG/`

### This file contains code and notes for
1) cleaning contaminants using tapioca
2) parsing barcodes
3) splitting fastqs 
4) 
6) 
7) 

## Cleaning contaminants

Being executed on ponderosa using tapioca pipeline. Commands in bash script, executed as below (10/16/21).

    $ module load fqutils/0.4.1
    $ module load bowtie2/2.2.5
    $ bash cleaning_bash.sh &

After GOAG.clean.fastq has been produced, clean out duplicate raw data:

    $ rm -rf GOAG-lib10_S1_L001_R1_001.fastq
 
Number of reads **before** cleaning:

    $ grep -c "^@" GOAG-lib10_S1_L001_R1_001.fastq > number_of_rawreads.txt
    $ less number_of_rawreads.txt
    # 2352951075
    
Number of reads **after** cleaning:

    $ grep -c "^@" GOAG.clean.fastq > number_of_cleanreads.txt
    $ less number_of_cleanreads.txt
    # 

## Barcode parsing:

Barcode keyfile is `/working/parchman/GOAG/XXXXXXXXX_bcode.csv`
  
    $ perl parse_barcodes768.pl barcodeKey_CLEAN_desertTortoises.csv GOAG.clean.fastq A00 &

`NOTE`: the A00 object is the code that identifies the sequencer (first three characters after the @ in the fastq identifier).

    $ less parsereport_GOAG.clean.fastq
    #Good mids count: 1549613635
    #Bad mids count: 52409522
    #Number of seqs with potential MSE adapter in seq: 448456
    #Seqs that were too short after removing MSE and beyond: 137


Cleaning up the directory:

    $ rm GOAG.clean.fastq
    $ rm miderrors_GOAG.clean.fastq
    $ rm parsereport_GOAG.clean.fastq


## Splitting fastq by individual ID

Make ids file

    $ cut -f 3 -d "," barcodeKey_CLEAN_desertTortoises.csv | grep "_" > GOAG_ids_noheader.txt
    # Note: 698 individuals


Split fastqs by individual, put in a new directory

    $ mkdir raw_fastqs
    $ perl splitFastq_universal_regex.pl GOAG_ids_noheader.txt parsed_GOAG.clean.fastq &

Zip the parsed*fastq files for now, but delete once patterns and qc are verified:

    $ gzip GOAG.clean.fastq

# Done to here 11/16/21

Total reads for GOAG (698 individuals)

    $ grep -c "^@" raw_fastqs/*fastq > seqs_per_ind.txt

Summarize in R

    R
    dat <- read.delim("seqs_per_ind.txt", header=F, sep=":")
        dim(dat)
        head(dat)
        
    sum(dat[,2])
        

Zip fastqs:

    $ gzip raw_fastqs/*fastq
