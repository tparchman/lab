#!/bin/bash

##########################################################################################
## MGS cleaning 
##########################################################################################

/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/illumina_oligos --pct 20 /working/parchman/MGS/MGS-1_S1_L001_R1_001.fastq > MGS.readstofilter.ill.txt 

echo "Illumina filtering done for lane 1"

/working/jahner/tapioca/src/tap_contam_analysis --db /archive/parchman_lab/rawdata_to_backup/contaminants/phix174 --pct 80 /working/parchman/MGS/MGS-1_S1_L001_R1_001.fastq > MGS.readstofilter.phix.txt 

echo "PhiX filtering done for lane 1"


/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/ecoli-k-12 --pct 80 /working/parchman/MGS/MGS-1_S1_L001_R1_001.fastq > MGS.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 1"


cat /working/parchman/MGS/MGS-1_S1_L001_R1_001.fastq | fqu_cull -r MGS.readstofilter.ill.txt MGS.readstofilter.phix.txt MGS.readstofilter.ecoli.txt > MGS.clean.fastq

echo "Clean copy of lane 1 done"


