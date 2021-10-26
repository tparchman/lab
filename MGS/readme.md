## Mojave Ground Squirrel analysis notes

### Notes on contaminant cleaning and barcode parsing 10-21

`NOTE`: One library was sequenced on one Novaseq lane in late September of 2021
`NOTE`: Contaminant cleaning and barcode parsing in `/working/parchman/MGS/`

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

After MGS.clean.fastq has been produced, clean out duplicate raw data:

    $ rm -rf MGS-1_S1_L001_R1_001.fastq
 
Number of reads **before** cleaning:

    $ grep -c "^@" MGS-1_S1_L001_R1_001.fastq > number_of_rawreads.txt &
    $ less number_of_rawreads.txt
    # 2,256,351,365
    
Number of reads **after** cleaning:

    $ grep -c "^@" MGS.clean.fastq > number_of_cleanreads.txt &
    $ less number_of_cleanreads.txt
    # 1,707,635,152
        
## Barcode parsing:

Barcode keyfile is `/working/parchman/MGS/barcodeKey_lib9_mojaveGroundSquirrels.csv`
  
    $ perl parse_barcodes768.pl barcodeKey_lib9_mojaveGroundSquirrels.csv MGS.clean.fastq A00 &

`NOTE`: the A00 object is the code that identifies the sequencer (first three characters after the @ in the fastq identifier).

    $ less parsereport_MGS.clean.fastq
    #Good mids count: 
    #Bad mids count: 
    #Number of seqs with potential MSE adapter in seq: 
    #Seqs that were too short after removing MSE and beyond: 
          
Cleaning up the directory:

    $ rm MGS.clean.fastq
    $ rm miderrors_MGS.clean.fastq
    $ rm parsereport_MGS.clean.fastq

## Splitting fastq by individual ID

Make ids file

    $ cut -f 3 -d "," barcodeKey_lib9_mojaveGroundSquirrels.csv | grep "_" > MGS_ids_noheader.txt
    # Note: 522 individuals

# Done to here 10.26.21

Split fastqs by individual, put in a new directory

    $ mkdir raw_fastqs
    $ perl splitFastq_universal_regex.pl MGS_ids_noheader.txt parsed_MGS.clean.fastq &

Zip the parsed*fastq files for now, but delete once patterns and qc are verified:

    $ gzip parsed_S1_11_20.clean.fastq


Total reads for mojave ground squirrels (XXX individuals)

    $ grep -c "^@" raw_fastqs/*fastq > seqs_per_ind.txt

Summarize in R

    R
    dat <- read.delim("seqs_per_ind.txt", header=F, sep=":")
        dim(dat)
        head(dat)
        
    sum(dat[,2])
        

Zip fastqs:

    $ gzip raw_fastqs/*fastq
