## Data available for analysis

Here I'm listing directories for a number of projects for which DNA sequencing has been completed. Trevor should be able to give you files with the geographic coordinates for each project.

1). Balsamorhiza sagittata (arrowleaf balsamroot)

Raw Novaseq data:

    /archive/parchman_lab/rawdata_to_backup/BASA_CHDO_rawNOVASEQ/Library2-BASA_S2_L002_R1_001.fastq.gz

Cleaned and demultiplexed .fastq files by individual:

    /working/tfaske/balsam/demult/fastq (BASA)
    /working/tfaske/balsam/demult/outgroup_fastq (outgroup)

Trevor's github page for BASA:
https://github.com/trevorfaske/BASA

2). Chaenactic douglassi (Douglas dusty maiden)

Raw Novaseq data:

    /archive/parchman_lab/rawdata_to_backup/BASA_CHDO_rawNOVASEQ/Library1-CHDO_S1_L001_R1_001.fastq.gz

Cleaned and demultiplexed .fastq files by individual:

    /working/tfaske/dusty/demult/fastq (CHDO)
    /working/tfaske/dusty/demult/outgroup_fastq 
    (outgroup)

3). Achnatherum thurberianum (Thurbers needlegrass)

Cleaned and demultiplexed .fastq files by individual:

    /archive/parchman_lab/rawdata_to_backup/GSAF_11_20_bc_parse/ACTH

4). Poa Secunda (sandberg bluegrass)

Cleaned and demultiplexed .fastq files by individual:

    /archive/parchman_lab/rawdata_to_backup/GSAF_11_20_bc_parse/POSE



## Analyses to read up on:
Genetic environment association (GEA) analyses (RDA, LFMM, Gradient Forests).

Estimating ploidy from GBS data (gbs2ploidy, from Zach Gompert, and any other methods that are out there). Ideally we use such a method with the C. douglasii data, for which we know ploidy (2, 4, and a few 6) for most populations, but not all.

Methods for inferring gene flow for landscape genomic data. We havent been using such methods, but may like to in the near future.