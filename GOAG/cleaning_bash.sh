#!/bin/bash

##########################################################################################
## GOAG cleaning 
##########################################################################################

/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/illumina_oligos --pct 20 /working/parchman/GOAG/GOAG-lib10_S1_L001_R1_001.fastq > GOAG.readstofilter.ill.txt 

echo "Illumina filtering done for lane 1"

/working/jahner/tapioca/src/tap_contam_analysis --db /archive/parchman_lab/rawdata_to_backup/contaminants/phix174 --pct 80 /working/parchman/GOAG/GOAG-lib10_S1_L001_R1_001.fastq > GOAG.readstofilter.phix.txt 

echo "PhiX filtering done for lane 1"


/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/ecoli-k-12 --pct 80 /working/parchman/GOAG/GOAG-lib10_S1_L001_R1_001.fastq > GOAG.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 1"


cat /working/parchman/GOAG/GOAG-lib10_S1_L001_R1_001.fastq | fqu_cull -r GOAG.readstofilter.ill.txt GOAG.readstofilter.phix.txt GOAG.readstofilter.ecoli.txt > GOAG.clean.fastq

echo "Clean copy of lane 1 done"


