# Welcome to the notes file. Let's crush it!

JPJ 7/31/20

`NOTE`: This dataset was sequenced on one Novaseq lane with a bunch of sagebrush samples. Will effectively be working on both datasets (N = 608) until I get to split fastqs.

`NOTE`: All analyses will be completed in /working/jahner/grouse2/

## What's in this notes file?
1) tapioca
2) parsing barcodes
3) splitting fastqs
4) steal genome from NCBI
5) reference-based assembly
6) calling variants
7) filtering
8) entropy
9) diversity


## 1) tapioca

[tapioca pipeline](https://github.com/ncgr/tapioca)

`INITIAL FILE LOCATION`: /archive/parchman_lab/rawdata_to_backup/7_20_Tpodura_Atridentata/JA20188-184712528/SA20102-L1_sagebrush-300180278/L1-sagebrush_S1_L001_R1_001.fastq.gz

    $ cp L1-sagebrush_S1_L001_R1_001.fastq.gz /working/jahner/grouse2/ &
    $ gunzip L1-sagebrush_S1_L001_R1_001.fastq.gz &

Number of reads before cleaning:

    $ grep -c "^@" L1-sagebrush_S1_L001_R1_001.fastq
            2136725522

`NOTE`: need to create and modify a new tapioca bash script for each project

    $ module load fqutils/0.4.1
    $ module load bowtie2/2.2.5
    $ bash cleaning_bash_sage.sh &

Number of reads after cleaning:

    $ grep -c "^@" sage.clean.fastq
            1615589404

Cleaning up the directory:

    $ rm L1-sagebrush_S1_L001_R1_001.fastq
    $ rm sage.readstofilter.*



## 2) parsing barcodes

Obtained unformatted barcode key file from Trevor (LibMar_1_ARTR_GRSG.csv). Fixed the really weird line endings in `TextWrangler` (LibMar_1_ARTR_GRSG_good_endings.csv). Pulled out the correct columns and changed header with `Unix`.

    $ cut -d ',' -f 3,4,5 LibMar_1_ARTR_GRSG_good_endings.csv | sed 's/BARCODE_ID,BARCODE,Name/BARCODE_ID,BARCODE,ID/g' > sage_barcode_key.csv

Parsing barcodes 

    $ perl /working/jahner/perl_scripts/parse_barcodes768.pl sage_barcode_key.csv sage.clean.fastq A00 &

`NOTE`: the A00 object is the code that identifies the sequencer (first three characters after the @ in the fastq identifier).

    $ less parsereport_sage.clean.fastq
            Good mids count: 1552129971
            Bad mids count: 63459273
            Number of seqs with potential MSE adapter in seq: 272045
            Seqs that were too short after removing MSE and beyond: 160

Cleaning up the directory:

    $ rm sage.clean.fastq
    $ rm miderrors_sage.clean.fastq
    $ rm parsereport_sage.clean.fastq


## 3) splitting fastqs

Make ids file

    $ cut -f 3 -d "," sage_barcode_key.csv | grep "_" > sage_ids_noheader.txt

Split fastqs by individual

    $ mkdir raw_fastqs
    $ cd raw_fastqs
    $ perl /working/jahner/perl_scripts/splitFastq_universal_regex.pl ../sage_ids_noheader.txt ../parsed_sage.clean.fastq &

`IMPORTANT NOTE`: one grouse individual `MA_NB_549` (plate 6, well F5) did not split at all, so it must have missed an ingredient during library prep. It's possible that the individuals above (E5: `MA_TL_14162`) or below (G5: `SH_LS_14323`) on the well might have gotten double adaptors, so double check they aren't goofy during analyses.

Clean up the directory:

    $ rm parsed_sage.clean.fastq

Total reads for grouse (127 individuals)

    $ grep -c "^@" raw_fastqs/*fastq > seqs_per_ind.txt

Summarize in R

    R
    dat <- read.delim("seqs_per_ind.txt", header=F, sep=":")
        dim(dat)
        head(dat)
        
    sum(dat[,2])
        236743153

Zip the fastqs

    $ gzip raw_fastqs/*fastq


## STOP!!!

All individuals fastqs need to be backed up on a lab hard drive before continuing. THIS IS VERY IMPORTANT



## 4) steal genome from NCBI

`NOTE`: will be using the Gunnison sage-grouse [assembly](https://www.ncbi.nlm.nih.gov/assembly/GCA_005890655.1/) as a reference (GCA_005890655.1)

Download genome from NCBI (can delete files other than the `fasta` if you want). Look up NCBI genome path [here](https://ftp.ncbi.nlm.nih.gov/genomes/)

    $ mkdir reference_genome
    $ cd reference_genome
    $ rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/005/890/655/GCA_005890655.1_Cmin_1.0/ .

Remove 60 line endings

    $ gunzip GCA_005890655.1_Cmin_1.0_genomic.fna.gz
    $ perl /working/jahner/perl_scripts/remove60_fasta.pl GCA_005890655.1_Cmin_1.0_genomic.fna

Change scaffold names to something shorter

    $ perl /working/jahner/perl_scripts/rename_scaff.pl no60_GCA_005890655.1_Cmin_1.0_genomic.fna
    $ mv renamed_no60_GCA_005890655.1_Cmin_1.0_genomic.fna.txt gunnison_genome.fna
    $ gzip GCA_005890655.1_Cmin_1.0_genomic.fna
    $ rm no60_GCA_005890655.1_Cmin_1.0_genomic.fna
    
    

## 5) reference-based assembly

`NOTE`: Moving forward with 127 individual grouse

Make index for `bwa`

    $ module load bwa/0.7.8
    $ bwa index -p gunnison -a bwtsw reference_genome/gunnison_genome.fna &

Perl wrapper to run `bwa`. `NOTE`: `runbwa.pl` needs to be modified for every project (index name, output directory, number of cpus). 

    $ module load bwa/0.7.8
    $ perl runbwa.pl raw_fastqs/*fastq &

Convert `sam` files to `bam` files. `NOTE`: change number of threads based on current server usage.

    $ module load samtools/1.3
    $ perl /working/jahner/perl_scripts/sam2bamV1.3.pl *.sam

Cleaning up the directory
    
    $ rm *.sam
    $ rm *.sai



## 6) calling variants

`NOTE`: the four lines below are one single command. You should only change the output name.

    $ module load samtools/1.3 
    $ module load bcftools/1.3
    
    $ samtools mpileup -P ILLUMINA --BCF --max-depth 100 --adjust-MQ 50 --min-BQ 20 --min-MQ 20 --skip-indels --output-tags AD,DP --fasta-ref /working/jahner/grouse2/reference_genome/gunnison_genome.fna aln*sorted.bam | \
    bcftools call -m --variants-only --format-fields GQ --skip-variants indels | \
    bcftools filter --set-GTs . -i 'QUAL > 19 && FMT/GQ >9' | \
    bcftools view -m 2 -M 2 -v snps --apply-filter "PASS" --output-type v --output-file variants_rawfiltered_4aug20.vcf

`NOTE`: once you are sure that the command is working, put the job in the background (`control z`; `bg`)



## 7) filtering

number of loci in vcf

    $ grep -c "^scaffold" variants_rawfiltered_4aug20.vcf
            250962

make id file for reheadering

    $ ls *sorted.bam > bams.txt
    $ sed -s "s/aln_//" bams.txt | sed -s "s/.sorted.bam//" > grouse2_ids_col.txt

reheader vcf

    $ module load bcftools/1.3
    $ module load vcftools/0.1.14
    $ bcftools reheader -s grouse2_ids_col.txt variants_rawfiltered_4aug20.vcf -o rehead_variants_rawfiltered_4aug20.vcf

initial round of filtering (just getting a feeling of how filtering parameters might shape the dataset)

    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf5_miss9 --remove-filtered-all --maf 0.05 --max-missing 0.9 --recode --thin 100
            After filtering, kept 2192 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf5_miss8 --remove-filtered-all --maf 0.05 --max-missing 0.8 --recode --thin 100
	        After filtering, kept 13895 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf5_miss7 --remove-filtered-all --maf 0.05 --max-missing 0.7 --recode --thin 100
	        After filtering, kept 23831 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf5_miss6 --remove-filtered-all --maf 0.05 --max-missing 0.6 --recode --thin 100
	        After filtering, kept 30933 out of a possible 250962 Sites (reference)

    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf4_miss9 --remove-filtered-all --maf 0.04 --max-missing 0.9 --recode --thin 100
	        After filtering, kept 2816 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf4_miss8 --remove-filtered-all --maf 0.04 --max-missing 0.8 --recode --thin 100
	        After filtering, kept 16527 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf4_miss7 --remove-filtered-all --maf 0.04 --max-missing 0.7 --recode --thin 100
	        After filtering, kept 27230 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf4_miss6 --remove-filtered-all --maf 0.04 --max-missing 0.6 --recode --thin 100
	        After filtering, kept 34642 out of a possible 250962 Sites (reference)

    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf3_miss9 --remove-filtered-all --maf 0.03 --max-missing 0.9 --recode --thin 100
	        After filtering, kept 3823 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf3_miss8 --remove-filtered-all --maf 0.03 --max-missing 0.8 --recode --thin 100
	        After filtering, kept 20111 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf3_miss7 --remove-filtered-all --maf 0.03 --max-missing 0.7 --recode --thin 100
	        After filtering, kept 31963 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf3_miss6 --remove-filtered-all --maf 0.03 --max-missing 0.6 --recode --thin 100
	        After filtering, kept 39795 out of a possible 250962 Sites (reference)

    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf2.5_miss9 --remove-filtered-all --maf 0.025 --max-missing 0.9 --recode --thin 100
	        After filtering, kept 4547 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf2.5_miss8 --remove-filtered-all --maf 0.025 --max-missing 0.8 --recode --thin 100
	        After filtering, kept 22157 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf2.5_miss7 --remove-filtered-all --maf 0.025 --max-missing 0.7 --recode --thin 100
	        After filtering, kept 34718 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf2.5_miss6 --remove-filtered-all --maf 0.025 --max-missing 0.6 --recode --thin 100
	        After filtering, kept 42779 out of a possible 250962 Sites (reference)

    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf2_miss9 --remove-filtered-all --maf 0.02 --max-missing 0.9 --recode --thin 100
	        After filtering, kept 5185 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf2_miss8 --remove-filtered-all --maf 0.02 --max-missing 0.8 --recode --thin 100
	        After filtering, kept 24306 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf2_miss7 --remove-filtered-all --maf 0.02 --max-missing 0.7 --recode --thin 100
	        After filtering, kept 37687 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf2_miss6 --remove-filtered-all --maf 0.02 --max-missing 0.6 --recode --thin 100
	        After filtering, kept 45974 out of a possible 250962 Sites (reference)

    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf1_miss9 --remove-filtered-all --maf 0.01 --max-missing 0.9 --recode --thin 100
	        After filtering, kept 6835 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf1_miss8 --remove-filtered-all --maf 0.01 --max-missing 0.8 --recode --thin 100
	        After filtering, kept 30074 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf1_miss7 --remove-filtered-all --maf 0.01 --max-missing 0.7 --recode --thin 100
	        After filtering, kept 45604 out of a possible 250962 Sites (reference)
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --out variants_maf1_miss6 --remove-filtered-all --maf 0.01 --max-missing 0.6 --recode --thin 100
	        After filtering, kept 54833 out of a possible 250962 Sites (reference)

`NOTE`: moving forward with `maf3_miss8`

filter out bad individuals, then refilter

    $ vcftools --vcf variants_maf3_miss8.recode.vcf --missing-indv

`OPTIONAL`: look at distribution of missing data in `R`

    R
    j <- read.delim("out.imiss", header=T)
        dim(j)
        head(j)
    hist(j[,5], breaks=100, col="gray")

`NOTE`: will remove individuals with >50% missing data (N = 12)

    $ mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
    $ vcftools --vcf rehead_variants_rawfiltered_4aug20.vcf --remove lowDP.indv --recode --recode-INFO-all --out rehead_variants_rawfiltered_4aug20_noBadInds
        After filtering, kept 115 out of 127 Individuals
        After filtering, kept 250962 out of a possible 250962 Sites

clean up steps

    $ rm variants_maf*
    $ gzip variants_rawfiltered_4aug20.vcf &
    $ gzip rehead_variants_rawfiltered_4aug20.vcf &

now filter `maf3_miss8` based on new set of individuals

    $ vcftools --vcf rehead_variants_rawfiltered_4aug20_noBadInds.recode.vcf --out variants_maf3_miss8 --remove-filtered-all --maf 0.03 --max-missing 0.8 --recode --thin 100
            After filtering, kept 29158 out of a possible 250962 Sites (reference)

now need to rekill bad inds (missing >50%) and see where I am at (3 inds removed)

    $ vcftools --vcf variants_maf3_miss8.recode.vcf --missing-indv
    $ mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
    $ vcftools --vcf variants_maf3_miss8.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out variants_maf3_miss8_noBadInds
			After filtering, kept 112 out of 115 Individuals
            After filtering, kept 29158 out of a possible 29158 Sites

generate mpgl

    $ perl /working/jahner/perl_scripts/vcf2mpglV1.3TLP.pl variants_maf3_miss8_noBadInds.recode.vcf

create ids file

    $ vcftools --vcf variants_maf3_miss8_noBadInds.recode.vcf --missing-indv 
    $ cut -f 1 out.imiss > grouse2_ids_112.txt
    $ sed "s/INDV/ind/" grouse2_ids_112.txt | sed "s/aln_//g" | sed "s/.sorted.bam//g" > grouse2_ids_112_good_head.txt

calculate coverage

    $ module load samtools/1.3
    $ module load bcftools/1.3
    $ perl /working/jahner/perl_scripts/coverage_calc_bam.pl variants_maf3_miss8_noBadInds.recode.mpgl grouse2_ids_112_good_head.txt /working/jahner/grouse2/sam_sai/ &
	
loc names

    $ cut -d " " -f 1 variants_maf3_miss8_noBadInds.recode.mpgl > loc_names_29158.txt

filter out over-assembled loci

    R
    dat <- read.csv("variants_maf3_miss8_noBadInds.recode.mpgl_coverage.csv", header=F)
		dim(dat)
		dat[1:10,1:10]

    dat_noname <- dat[,-1]
		dim(dat_noname)
		dat_noname[1:10,1:10]

    loc_names <- read.delim("loc_names_29158.txt", header=F)
		head(loc_names)

    avg_29158 <- vector()
    in_out_29158_15 <- vector()
    in_out_29158_20 <- vector()
    in_out_29158_25 <- vector()
    in_out_29158_30 <- vector()

    for (i in 1:29158)
		{
		avg <- mean(dat_noname[,i])
		avg_29158 <- append(avg_29158, avg)
	
		if (avg <= 15)	{ in_out_29158_15 <- append(in_out_29158_15, 1) }
		else			{ in_out_29158_15 <- append(in_out_29158_15, 0) }

		if (avg <= 20)	{ in_out_29158_20 <- append(in_out_29158_20, 1) }
		else			{ in_out_29158_20 <- append(in_out_29158_20, 0) }
	
		if (avg <= 25)	{ in_out_29158_25 <- append(in_out_29158_25, 1) }
		else			{ in_out_29158_25 <- append(in_out_29158_25, 0) }
	
		if (avg <= 30)	{ in_out_29158_30 <- append(in_out_29158_30, 1) }
		else			{ in_out_29158_30 <- append(in_out_29158_30, 0) }
	
		}

    sum(in_out_29158_15) / 29158
		26566 (92.0%) - (reference)
    sum(in_out_29158_20) / 29158
		28261 (96.9%) - (reference)
    sum(in_out_29158_25) / 29158
		28755 (97.8%) - (reference)
    sum(in_out_29158_30) / 29158
		28946 (98.4%) - (reference)

    sub_20 <- dat_noname[,in_out_29158_20==1]
		dim(sub_20)
    sub_20_avg <- subset(avg_29158, in_out_29158_20==1)	

    kill_locs <- subset(loc_names, in_out_29158_20==0)
    dim(kill_locs)
		head(kill_locs)

    write.table(kill_locs, file="high_cov_loc_list_to_be_removed.txt", quote=F, row.names=F, col.names=F)

actually filter in unix. `NOTE`: removing loci with mean coverage > 20

    $ sed "s/:/\t/" high_cov_loc_list_to_be_removed.txt > high_cov_loc_list_to_be_removed_tabdelim.txt
	$ module load vcftools/0.1.14 
	$ vcftools --vcf variants_maf3_miss8_noBadInds.recode.vcf --exclude-positions high_cov_loc_list_to_be_removed_tabdelim.txt --recode --recode-INFO-all --out variants_maf3_miss8_noBadInds_noHighCov
			After filtering, kept 28261 out of a possible 29158 Sites

final filter on FIS

	$ gzip variants_maf3_miss8_noBadInds_noHighCov.recode.vcf
	$ ~/miniconda3/bin/python3.6 /working/jahner/python_scripts/Fis_filtering_vcf.py /working/jahner/grouse2/sam_sai/ variants_maf3_miss8_noBadInds_noHighCov.recode.vcf.gz
			TOTAL LOCI:	 28261
			loci Fis < -0.5: 56
			loci Fis > 0.5: 197
			Filtered loci:	 56
			REMAINING LOCI:	 28205
	$ rm Fis_filtered.recode.vcf.gz 

calculate final coverage

	$ module load samtools/1.3
	$ module load bcftools/1.3
	$ perl /working/jahner/perl_scripts/vcf2mpglV1.3TLP.pl Fis_filtered.recode.vcf
	$ perl /working/jahner/perl_scripts/coverage_calc_bam.pl Fis_filtered.recode.mpgl grouse2_ids_112_good_head.txt /working/jahner/grouse2/sam_sai/ &

	R
	j <- read.csv("Fis_filtered.recode.mpgl_coverage.csv", header=F)
	k <- j[,-1]
	mean_vect <- vector()
	for (i in 1:112) { mean_vect <- append(mean_vect, mean(as.numeric(k[i,]))) }
	mean(mean_vect)
			8.318152

cleaning up the directory

	$ gzip aln* &
	$ rm variants_maf3_miss8_noBadInds_noHighCov.recode.vcf.gz.*
	$ gzip *coverage.csv
	$ gzip *.vcf &
	$ gzip variants_maf3_miss8_noBadInds.recode.mpgl



## 8) entropy

generate pntest file

	$ perl /working/jahner/perl_scripts/gl2genestV1.3.pl Fis_filtered.recode.mpgl mean

generate pops file

	$ cut -d "_" -f 1,2 grouse2_ids_112_good_head.txt > grouse2_pops_112.txt

make genotype likelihood matrix

	R
	read.table("pntest_mean_Fis_filtered.recode.txt", header=F)->gl
	read.table("grouse2_ids_112_good_head.txt", header=T)->ids
	read.table("grouse2_pops_112.txt", header=T)->pops
	t(gl)->tgl
	cbind(ids, pops, tgl)->tidsgl
	write.table(tidsgl, file="grouse_gl_matrix_112_28205.txt", sep=" ", row.names=F, col.names=F , quote=F)

generate lda files

	R
	g <- read.table("pntest_mean_Fis_filtered.recode.txt", header=F)
	names <- read.table("grouse2_ids_112.txt", header=T)
	pops <- read.table("grouse2_pops_112.txt", header=T)
	nind <- dim(g)[2]
	nloci <- dim(g)[1]

	gmn<-apply(g,1,mean, na.rm=T)
	gmnmat<-matrix(gmn,nrow=nloci,ncol=nind)
	gprime<-g-gmnmat ## remove mean
	gcovarmat<-matrix(NA,nrow=nind,ncol=nind)
	for(i in 1:nind){
	    for(j in i:nind){
        	if (i==j){
            	gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
        	}
        	else{
            	gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
            	gcovarmat[j,i]<-gcovarmat[i,j]
        	}
    	}
	}

	pcgcov<-prcomp(x=gcovarmat,center=TRUE,scale=FALSE)
	pcgcov->pcg
	library(MASS)
	k2<-kmeans(pcg$x[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
	k3<-kmeans(pcg$x[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
	k4<-kmeans(pcg$x[,1:5],4,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
	k5<-kmeans(pcg$x[,1:5],5,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
	k6<-kmeans(pcg$x[,1:5],6,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
	k7<-kmeans(pcg$x[,1:5],7,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

	ldak2<-lda(x=pcg$x[,1:5],grouping=k2$cluster,CV=TRUE)
	ldak3<-lda(x=pcg$x[,1:5],grouping=k3$cluster,CV=TRUE)
	ldak4<-lda(x=pcg$x[,1:5],grouping=k4$cluster,CV=TRUE)
	ldak5<-lda(x=pcg$x[,1:5],grouping=k5$cluster,CV=TRUE)
	ldak6<-lda(x=pcg$x[,1:5],grouping=k6$cluster,CV=TRUE)
	ldak7<-lda(x=pcg$x[,1:5],grouping=k7$cluster,CV=TRUE)

	write.table(round(ldak2$posterior,5),file="ldak2.txt",quote=F,row.names=F,col.names=F)
	write.table(round(ldak3$posterior,5),file="ldak3.txt",quote=F,row.names=F,col.names=F)
	write.table(round(ldak4$posterior,5),file="ldak4.txt",quote=F,row.names=F,col.names=F)
	write.table(round(ldak5$posterior,5),file="ldak5.txt",quote=F,row.names=F,col.names=F)
	write.table(round(ldak6$posterior,5),file="ldak6.txt",quote=F,row.names=F,col.names=F)
	write.table(round(ldak7$posterior,5),file="ldak7.txt",quote=F,row.names=F,col.names=F)

make mpgl input file for entropy

	$ grep "_" grouse2_ids_112_good_head.txt > grouse2_ids_112_nohead.txt
	$ perl /working/jahner/perl_scripts/create_entropy_top_2rows.pl grouse2_ids_112_nohead.txt 
	$ mkdir entropy
	$ mv ldak* entropy/
	$ mv entropy_2rows.txt entropy/
	$ cd entropy/
	$ cat entropy_2rows.txt ../Fis_filtered.recode.mpgl > grouse2_entropy.mpgl
	
`NOTE`: need to add a header to the top of the mpgl with nano: `112 28205 1`

running entropy

	$ module load entropy/1.2
	$ entropy -i grouse2_entropy.mpgl -o grouse2_k2_rep4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 2 -q ldak2.txt -m 1 -w 0 &> k2stdout_rep4.txt &
	$ entropy -i grouse2_entropy.mpgl -o grouse2_k3_rep4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 3 -q ldak3.txt -m 1 -w 0 &> k3stdout_rep4.txt &
	$ entropy -i grouse2_entropy.mpgl -o grouse2_k4_rep4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 4 -q ldak4.txt -m 1 -w 0 &> k4stdout_rep4.txt &
	$ entropy -i grouse2_entropy.mpgl -o grouse2_k5_rep4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 5 -q ldak5.txt -m 1 -w 0 &> k5stdout_rep4.txt &
	$ entropy -i grouse2_entropy.mpgl -o grouse2_k6_rep4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 6 -q ldak6.txt -m 1 -w 0 &> k6stdout_rep4.txt &
	$ entropy -i grouse2_entropy.mpgl -o grouse2_k7_rep4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak7.txt -m 1 -w 0 &> k7stdout_rep4.txt &

entropy DICs

	$ module load entropy/1.2
	$ estpost.entropy grouse2_k2_rep4.hdf5 -s 3 -p deviance
	$ estpost.entropy grouse2_k3_rep4.hdf5 -s 3 -p deviance
	$ estpost.entropy grouse2_k4_rep4.hdf5 -s 3 -p deviance
	$ estpost.entropy grouse2_k5_rep4.hdf5 -s 3 -p deviance
	$ estpost.entropy grouse2_k6_rep4.hdf5 -s 3 -p deviance
	$ estpost.entropy grouse2_k7_rep4.hdf5 -s 3 -p deviance

generate entropy q files

	$ module load entropy/1.2
	$ estpost.entropy grouse2_k2_rep3.hdf5  -p q -s 0 -o q2_rep3.txt
	$ estpost.entropy grouse2_k3_rep3.hdf5  -p q -s 0 -o q3_rep3.txt
	$ estpost.entropy grouse2_k4_rep3.hdf5  -p q -s 0 -o q4_rep3.txt
	$ estpost.entropy grouse2_k5_rep3.hdf5  -p q -s 0 -o q5_rep3.txt
	$ estpost.entropy grouse2_k6_rep3.hdf5  -p q -s 0 -o q6_rep3.txt

generate gprobs file

	$ module load entropy/1.2
	$ estpost.entropy  grouse2_k2_rep0.hdf5 -p gprob -s 0 -o gprob2_rep0.txt



## 9) diversity

make bam lists per pop

	$ grep "HA_HA" grouse2_ids_112.txt > HA_HA_bams.txt
	$ grep "MA_FA" grouse2_ids_112.txt > MA_FA_bams.txt
	$ grep "MA_NB" grouse2_ids_112.txt > MA_NB_bams.txt
	$ grep "MA_TL" grouse2_ids_112.txt > MA_TL_bams.txt
	$ grep "SH_HL" grouse2_ids_112.txt > SH_HL_bams.txt
	$ grep "SH_LS" grouse2_ids_112.txt > SH_LS_bams.txt


	$ sed "s/\$/\.sorted\.bam/g" HA_HA_bams.txt | sed "s/^/aln_/g" > HA_HA_bam_names.txt
	$ sed "s/\$/\.sorted\.bam/g" MA_FA_bams.txt | sed "s/^/aln_/g" > MA_FA_bam_names.txt
	$ sed "s/\$/\.sorted\.bam/g" MA_NB_bams.txt | sed "s/^/aln_/g" > MA_NB_bam_names.txt
	$ sed "s/\$/\.sorted\.bam/g" MA_TL_bams.txt | sed "s/^/aln_/g" > MA_TL_bam_names.txt
	$ sed "s/\$/\.sorted\.bam/g" SH_HL_bams.txt | sed "s/^/aln_/g" > SH_HL_bam_names.txt
	$ sed "s/\$/\.sorted\.bam/g" SH_LS_bams.txt | sed "s/^/aln_/g" > SH_LS_bam_names.txt

first step (doSAF)

Options:

-bam <INPUT> = input of bam names for population
	
-doSaf <INPUT> = option 1: calculate the site allele frequency likelihood based on individual genotype likelihoods assuming HWE
	
-anc <INPUT> = ancestral fasta file (i.e. the genome)
	
-GL <INPUT> = genotype likelihood model (1 = SAMtools; 2 = GATK; 3 = SOAPsnp; 4 = SYK)
	
-P <INPUT> = number of cores used
	
-out <INPUT> = outfile prefix


	$ source activate angsd
	$ angsd -bam HA_HA_bam_names.txt -doSaf 1 -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 -P 2 -out HA_HA &
	$ angsd -bam MA_FA_bam_names.txt -doSaf 1 -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 -P 2 -out MA_FA &
	$ angsd -bam MA_NB_bam_names.txt -doSaf 1 -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 -P 2 -out MA_NB &
	$ angsd -bam MA_TL_bam_names.txt -doSaf 1 -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 -P 2 -out MA_TL &
	$ angsd -bam SH_HL_bam_names.txt -doSaf 1 -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 -P 2 -out SH_HL &
	$ angsd -bam SH_LS_bam_names.txt -doSaf 1 -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 -P 2 -out SH_LS &


Second step (realSFS): estimates the (multi) SFS based on a .saf.idx file generated from step 1 (doSaf)

-P <INPUT> = number of cores used

-fold = specify folded SFS (1)

	$ source activate angsd
	$ realSFS HA_HA.saf.idx -P 2 -fold 1 > HA_HA.sfs &
	$ realSFS MA_FA.saf.idx -P 2 -fold 1 > MA_FA.sfs &
	$ realSFS MA_NB.saf.idx -P 2 -fold 1 > MA_NB.sfs &
	$ realSFS MA_TL.saf.idx -P 2 -fold 1 > MA_TL.sfs &
	$ realSFS SH_HL.saf.idx -P 2 -fold 1 > SH_HL.sfs &
	$ realSFS SH_LS.saf.idx -P 2 -fold 1 > SH_LS.sfs &



Third step (doThetas): calculate the thetas (population scale mutation rates) for each site

-bam <INPUT> = input of bam names for population
	
-out <INPUT> = outfile prefix
	
-doThetas 1 = calculate the thetas (use 1, not sure what other option values exist or do)

-doSaf <INPUT> = option 1: calculate the site allele frequency likelihood based on individual genotype likelihoods assuming HWE
	
-pest = if -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allele frequency for each site

-anc <INPUT> = ancestral fasta file (i.e. the genome)
	
-GL <INPUT> = genotype likelihood model (1 = SAMtools; 2 = GATK; 3 = SOAPsnp; 4 = SYK)

	$ source activate angsd
	$ angsd -bam HA_HA_bam_names.txt -out HA_HA -doThetas 1 -doSaf 1 -pest HA_HA.sfs -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 & 
	$ angsd -bam MA_FA_bam_names.txt -out MA_FA -doThetas 1 -doSaf 1 -pest MA_FA.sfs -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 & 
	$ angsd -bam MA_NB_bam_names.txt -out MA_NB -doThetas 1 -doSaf 1 -pest MA_NB.sfs -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 & 
	$ angsd -bam MA_TL_bam_names.txt -out MA_TL -doThetas 1 -doSaf 1 -pest MA_TL.sfs -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 & 
	$ angsd -bam SH_HL_bam_names.txt -out SH_HL -doThetas 1 -doSaf 1 -pest SH_HL.sfs -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 & 
	$ angsd -bam SH_LS_bam_names.txt -out SH_LS -doThetas 1 -doSaf 1 -pest SH_LS.sfs -anc /working/jahner/grouse2/reference_genome/gunnison_genome.fna -GL 1 & 


Fourth step (thetaStat do_stat): estimate Tajimas D and other statistics

-win = window size

-step = step size (step must be equal or greater than win if you want to avoid overlap)

`NOTE`: D stored in *thetaWindow.pestPG

	$ source activate angsd
	$ thetaStat do_stat HA_HA.thetas.idx -win 50000 -step 50000 -outnames HA_HA.thetaWindow
	$ thetaStat do_stat MA_FA.thetas.idx -win 50000 -step 50000 -outnames MA_FA.thetaWindow
	$ thetaStat do_stat MA_NB.thetas.idx -win 50000 -step 50000 -outnames MA_NB.thetaWindow
	$ thetaStat do_stat MA_TL.thetas.idx -win 50000 -step 50000 -outnames MA_TL.thetaWindow
	$ thetaStat do_stat SH_HL.thetas.idx -win 50000 -step 50000 -outnames SH_HL.thetaWindow
	$ thetaStat do_stat SH_LS.thetas.idx -win 50000 -step 50000 -outnames SH_LS.thetaWindow


Fourth step data extraction in R (Tajima's D)

	R

	library(data.table)

	CI <- function(data){
	    low <- mean(data) - 1.96*(sd(data)/sqrt(length(data)))
	    high <- mean(data) + 1.96*(sd(data)/sqrt(length(data)))
 	   return(c(low,high))
 	}

	pestPG_files <- list.files(pattern='*Window.pestPG')
  	TajD_df <- as.data.frame(matrix(nrow=length(pestPG_files),ncol=4))
	names(TajD_df) <- c('Pop','TajD','TajD_low','TajD_high')
  
	for (i in 1:length(pestPG_files)){
    	infile <- fread(pestPG_files[i])
    	Pop <- unlist(strsplit(pestPG_files[i],'.',fixed=TRUE))[1]
    	TajD_df$Pop[i] <- Pop
    	TajD <- infile$Tajima
    	TajD_df$TajD[i] <- mean(TajD)  
    	TajD_ci <- CI(TajD)
    	TajD_df$TajD_low[i] <- TajD_ci[1]
    	TajD_df$TajD_high[i] <- TajD_ci[2]
    }
	rm('infile')
	write.csv(TajD_df,'angsd_d_out.csv',row.names=FALSE)


Fifth step (thetaStat print)

`OUTPUTS`: chromosome; position; Theta Watterson; Nucleotide diversity; Theta Singleton category; Theta H; Theta L

	$ source activate angsd
	$ thetaStat print HA_HA.thetas.idx > HA_HA.theta_out
	$ thetaStat print MA_FA.thetas.idx > MA_FA.theta_out
	$ thetaStat print MA_NB.thetas.idx > MA_NB.theta_out
	$ thetaStat print MA_TL.thetas.idx > MA_TL.theta_out
	$ thetaStat print SH_HL.thetas.idx > SH_HL.theta_out
	$ thetaStat print SH_LS.thetas.idx > SH_LS.theta_out


Fifth step data extraction in R (pi and theta W)

	R

	library(data.table)

	CI <- function(data){
	  low <- mean(data) - 1.96*(sd(data)/sqrt(length(data)))
	  high <- mean(data) + 1.96*(sd(data)/sqrt(length(data)))
	  return(c(low,high))
	}

	subsetTheta_files <- list.files(pattern='*theta_out')
	subsetPisum_df <- as.data.frame(matrix(nrow=length(subsetTheta_files),ncol=7))
	names(subsetPisum_df) <- c('Pop','Pi','Pi_low','Pi_high', 'Watt', 'Watt_low', 'Watt_high')
	subsetPisamp_df <- as.data.frame(matrix(ncol=2))
	names(subsetPisamp_df) <- c('Pop','Pi')

	for (i in 1:length(subsetTheta_files)){
  		infile <- fread(subsetTheta_files[i])
  		Pop <- unlist(strsplit(subsetTheta_files[i],'.',fixed=TRUE))[1]
  		pi <- exp(infile$Pairwise) #in log form so must exp
  		watt <- exp(infile$Watterson)
  
  		subsetPisum_df$Pop[i] <- Pop
  		subsetPisum_df$Pi[i] <- mean(pi) #in log form 
  		subsetPisum_df$Watt[i] <- mean(watt) #in log form 
  
  		ci_pi <- CI(pi)
  		subsetPisum_df$Pi_low[i] <- ci_pi[1]
  		subsetPisum_df$Pi_high[i] <- ci_pi[2]
  
  		ci_watt <- CI(watt)
  		subsetPisum_df$Watt_low[i] <- ci_watt[1]
  		subsetPisum_df$Watt_high[i] <- ci_watt[2]
	}
	rm('infile')
	write.csv(subsetPisum_df,'angsd_piWatt_out.csv',row.names=FALSE)





