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
    # 
    
Number of reads **after** cleaning:

    $ grep -c "^@" GOAG.clean.fastq > number_of_cleanreads.txt
    $ less number_of_cleanreads.txt
    # 

# Done to here 11/9/21

## Barcode parsing:

Barcode keyfile is `/working/parchman/GOAG/XXXXXXXXX_bcode.csv`
  
    $ perl parse_barcodes768.pl final_timema3_Piper1_bcode.csv MGS.clean.fastq A00 &

`NOTE`: the A00 object is the code that identifies the sequencer (first three characters after the @ in the fastq identifier).

    $ less parsereport_tpodura.clean.fastq
    #Good mids count: 1617562664
    #Bad mids count: 58132389
    #Number of seqs with potential MSE adapter in seq: 305112
    #Seqs that were too short after removing MSE and beyond: 193
          
Cleaning up the directory:

    $ rm tpodura.clean.fastq
    $ rm miderrors_tpodura.clean.fastq
    $ rm parsereport_tpodura.clean.fastq
    
Total reads for T. podura (598 individuals)

    $ grep -c "^@" raw_fastqs/*fastq > seqs_per_ind.txt

Summarize in R

    R
    dat <- read.delim("seqs_per_ind.txt", header=F, sep=":")
        dim(dat)
        head(dat)
        
    sum(dat[,2])
        

Zip fastqs:

    $ gzip raw_fastqs/*fastq
