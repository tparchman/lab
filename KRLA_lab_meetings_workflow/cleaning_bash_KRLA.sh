#!/bin/bash

##########################################################################################
## FRLA2 cleaning 
##########################################################################################

/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/illumina_oligos --pct 20 /working/parchman/KRLA/KRLA_S1_L001_R1_001.fastq > KRLA.readstofilter.ill.txt 

echo "Illumina filtering done for lane 1"

/working/jahner/tapioca/src/tap_contam_analysis --db /archive/parchman_lab/rawdata_to_backup/contaminants/phix174 --pct 80 /working/parchman/KRLA/KRLA_S1_L001_R1_001.fastq > KRLA.readstofilter.phix.txt 

echo "PhiX filtering done for lane 1"


/working/jahner/tapioca/src/tap_contam_analysis --db  /archive/parchman_lab/rawdata_to_backup/contaminants/ecoli-k-12 --pct 80 /working/parchman/KRLA/KRLA_S1_L001_R1_001.fastq > KRLA.readstofilter.ecoli.txt

echo "ecoli filtering done for lane 1"


cat /working/parchman/KRLA/KRLA_S1_L001_R1_001.fastq | fqu_cull -r KRLA.readstofilter.ill.txt KRLA.readstofilter.phix.txt KRLA.readstofilter.ecoli.txt > KRLA.clean.fastq

echo "Clean copy of lane 1 done"

