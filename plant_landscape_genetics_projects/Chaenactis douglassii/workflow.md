### SNP calling
##### Generating the reference with dDocent:


```{r eval=FALSE}
dDocent:

conda activate ddocent_env

1. Create a new folder called RefOpt, with a subset of individuals of your total data set (less than 5)

./ReferenceOpt.sh 4 10 4 10 SE 5 0.95 0.99 0.005

2. Visualize data in kopt.data: plot values for each k1,k2 combination across similarity thresholds, pick a similarity threshold at the point of inflection on the curve 
.
****
library(ggplot2)
data.table <- read.table("kopt.data", header = FALSE, col.names= c("k1","k2","Similarity", "Contigs"))
data.table$K1K2 <- paste(data.table$k1, data.table$k2, sep=",")
df=data.frame(data.table)
df$K1K2 <- as.factor(df$K1K2)
p <- ggplot(df, aes(x=Similarity, y=Contigs, group=K1K2)) + scale_x_continuous(breaks=seq(0.8,0.98,0.01)) + geom_line(aes(colour = K1K2))
p
ggsave("kvalues.pdf",p,height=8,width = 10,units = 'in')
*****
---0.975, k 7 k 6---

3. Trimming the reads 

4. Run RefMapOpt.sh using the similarity threshold picked from last step 

./RefMapOpt.sh 6 7 6 7 0.975 SE 20 

5. Copy the reference to the main directory 

6. Trimming all the reads

7. Run dDocent on your full data set for reference mapping  97.5 7 6 
./RefMapOpt.sh 6 7 6 7 0.975 SE 20 
```

##### Calling SNPs with bcftools:

```{r eval=FALSE}
./bwa_sort.sh

*****
#!/bin/sh
INDS=($(for i in /working/mascaro/chdo/final/main/*.F.fq.gz; do echo $(basename ${i%.F.fq.gz*}); done))
for IND in ${INDS[@]};
do
	# declare variables
	REF=/working/mascaro/chdo/final/main/reference.fasta
	FORWARD=/working/mascaro/chdo/final/main/${IND}.F.fq.gz
	OUTPUT=/working/mascaro/chdo/final/main/${IND}_sort.bam
	# then align and sort
	echo "Aligning $IND with bwa"
	bwa mem -M -t 10 $REF $FORWARD \
	 | samtools view -b | \
	samtools sort -T ${IND} > $OUTPUT
done
*****

./bcftools.sh ****

*****
#!/bin/sh
REF=/working/mascaro/chdo/final/main/reference.fasta
bcftools mpileup -a AD,DP,SP -Ou -f reference.fasta ./*_sort.bam | bcftools call -f GQ,GP -mO z -o ./chdo.vcf.gz
done
```

##### Filtering:


```{r eval=FALSE}
1. Keep only Biallelic:

vcftools --remove-indels --min-alleles 2 --max-alleles 2 --remove-filtered-all  --recode --recode-INFO-all --gzvcf chdo.vcf.gz --out chdo_biallelic

#### After filtering, kept 276 out of 276 Individuals
#### Outputting VCF file...
#### After filtering, kept 665153 out of a possible 3082476 Sites

***
2. Remove by MAF, missing, and Thin:

vcftools --max-missing 0.6 --maf 0.02 --thin 100 --remove-filtered-all --recode --recode-INFO-all --gzvcf chdo_biallelic.recode.vcf --out chdo_miss60_thin100_maf2
###After filtering, kept 276 out of 276 Individuals
###Outputting VCF file...
###After filtering, kept 28471 out of a possible 665153 Sites
vcftools --max-missing 0.7 --maf 0.02 --thin 100 --remove-filtered-all --recode --recode-INFO-all --gzvcf chdo_biallelic.recode.vcf --out chdo_miss70_thin100_maf2
### After filtering, kept 276 out of 276 Individuals
### Outputting VCF file...
###After filtering, kept 26457 out of a possible 665153 Sites
vcftools --max-missing 0.7 --maf 0.03 --thin 100 --remove-filtered-all --recode --recode-INFO-all --gzvcf chdo_biallelic.recode.vcf --out chdo_miss70_thin100_maf3
###After filtering, kept 276 out of 276 Individuals
###Outputting VCF file...
###After filtering, kept 25628 out of a possible 665153 Sites
vcftools --max-missing 0.8 --maf 0.03 --thin 100 --remove-filtered-all --recode --recode-INFO-all --gzvcf chdo_biallelic.recode.vcf --out chdo_miss80_thin100_maf3
###After filtering, kept 276 out of 276 Individuals
###Outputting VCF file...
###After filtering, kept 22788 out of a possible 665153 Sites
vcftools --max-missing 0.8 --maf 0.04 --thin 100 --remove-filtered-all --recode --recode-INFO-all --gzvcf chdo_biallelic.recode.vcf --out chdo_miss80_thin100_maf3
###After filtering, kept 276 out of 276 Individuals
###Outputting VCF file...
###After filtering, kept 22037 out of a possible 665153 Sites

3. Calculate missing data using vcfR:

### Missing data
require(readr)
library(data.table)
require(MASS)
require(ggplot2)
library(vcfR)

vcf <- read.vcfR("chdo_final.vcf", verbose = FALSE)

#get positions
chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
pos_ID <- paste(chrom,pos,sep = ':')

#get pl 
dp <- extract.gt(vcf, element = 'DP')

## check out PL and pos_ID
print(length(pos_ID))
print(pos_ID[1:10])
dp[1:5,1:5]
str(dp[1:5,1:5])

#Calculate missing loci and individual
nloci <- ncol(dp)
nindv <- nrow(dp)
print(nloci)
print(nindv)

miss_loci <- apply(dp,1, function(d) length(which(d == "0"))/nloci)
print(summary(miss_loci))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.003623 0.028986 0.055332 0.097826 0.199275 

print(length(which(miss_loci <= .10))) #16647
print(length(which(miss_loci <= .20))) #22037
print(length(which(miss_loci <= .30))) #22037
print(length(which(miss_loci <= .40))) #22037
print(length(which(miss_loci <= .50))) #22037

#keep_miss <- pos_ID[which(miss_loci <= .20)]
#length(keep_miss)
#
#keep_miss20_df <- data.frame(chrom=sapply(keep_miss,function(s) unlist(strsplit(as.character(s),':'))[1]),
#                             pos=sapply(keep_miss,function(s) unlist(strsplit(as.character(s),':'))[2]))
#print(head(keep_miss20_df))    
#write.table(keep_miss20_df,'keep_miss20.txt',quote=F,row.names=F,col.names=F,sep='\t')
#vcftools --vcf chdo_final.vcf --remove-filtered-all --recode --recode-INFO-all --positions "keep_miss20.txt" --out chdo_final.vcf

4. Remove bad_indv with mean depth lower than 5

#CD_SA_5
#CD_SO_1

vcftools --gzvcf  chdo_ultimo.vcf.recode.vcf --remove-indels  --remove-filtered-all --recode --recode-INFO-all --remove bad_indvi.txt --out chdo_very_last
#After filtering, kept 274 out of 276 Individuals
#Outputting VCF file...
#After filtering, kept 22037 out of a possible 22037 Sites

mv chdo_very_last.recode.vcf chdo_ULTIMO.vcf


5. vcf statistics:
  
vcftools --gzvcf chdo_ULTIMO.vcf --site-quality --out quality
vcftools --gzvcf chdo_ULTIMO.vcf --freq2 --out  --max-alleles 2
vcftools --gzvcf chdo_ULTIMO.vcf --depth --out meandepthind
vcftools --gzvcf chdo_ULTIMO.vcf --site-mean-depth --out meandepsite
vcftools --gzvcf chdo_ULTIMO.vcf --missing-indv --out missing
vcftools --gzvcf chdo_ULTIMO.vcf --missing-site --out missingsite
vcftools --gzvcf chdo_ULTIMO.vcf --het --out het
vcftools --gzvcf chdo_ULTIMO.vcf --counts --out count


6. Recalculate missing data:

vcf <- read.vcfR("chdo_ULTIMO.vcf", verbose = FALSE)

#get positions
chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
pos_ID <- paste(chrom,pos,sep = ':')

dp <- extract.gt(vcf, element = 'DP')

print(length(pos_ID))
print(pos_ID[1:10])
dp[1:5,1:5]

ID <- colnames(dp)
ssp_ploid <- as.character(sapply(ID,function(s) unlist(strsplit(as.character(s),'_'))[1]))
ploidy <- sapply(ssp_ploid,function(s) gsub('(\\D)','',s,perl=TRUE))
print(length(ID))
print(ID[1:10])
print(ploidy[1:10])
print(unique(ploidy))

## Calculate missing loci and individual

print(dim(dp))
nindv <- ncol(dp)
nloci <- nrow(dp)
print(nloci)
print(nindv)

miss_loci <- apply(dp,1, function(d) length(which(d == "0"))/nindv)
#miss_loci <- apply(dp[,1:100],1, function(d) (length(which(d == "0"))/100))
print(length(miss_loci))
print(summary(miss_loci))

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.02555 0.05231 0.09489 0.20073 

###calculaate missing within each ploidy

#subset dp
dp2 <- dp[,which(ploidy=='2')]
dp4 <- dp[,which(ploidy=='4')]
dp6 <- dp[,which(ploidy=='6')]


#calc miss
miss_loci_2 <- apply(dp2,1, function(d) length(which(d == "0"))/ncol(dp2))
print(summary(miss_loci_2))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.04902 0.08530 0.14706 0.47059 

miss_loci_4 <- apply(dp4,1, function(d) length(which(d == "0"))/ncol(dp4))
print(summary(miss_loci_4))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.000000 0.007042 0.034429 0.056338 0.28873

miss_loci_6 <- apply(dp6,1, function(d) length(which(d == "0"))/ncol(dp6))
print(summary(miss_loci_6))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00000 0.02477 0.03333 0.5333 


miss_index_2 <- which(miss_loci_2 >= 0.95)
miss_index_4 <- which(miss_loci_4 >= 0.95)
miss_index_6 <- which(miss_loci_6 >= 0.95)

miss_index <- c(miss_index_2,miss_index_4,miss_index_6)
print(miss_index)

keep_missP <- pos_ID[-miss_index]
print(length(pos_ID))
print(length(miss_index))
print(length(keep_missP))


keep_missP_df <- data.frame(chrom=sapply(keep_missP,function(s) unlist(strsplit(as.character(s),':'))[1]),
                            pos=sapply(keep_missP,function(s) unlist(strsplit(as.character(s),':'))[2]))
print(head(keep_missP_df))    
# data frame with 0 columns and 0 rows

## Missing Inv
nloci <- nrow(dp)
#print(nloci)

miss_indv <- apply(dp,2, function(d) length(which(d == "0"))/nloci)
print(length(miss_indv))
print(summary(miss_indv))

#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003176 0.017947 0.039683 0.052309 0.078584 0.204701 

for (perc in seq(.1,.9,by=.1)){
  print(paste0('number of individuals with greater than ',perc*100,'% missing data.....  ',
               length(which(miss_indv >= perc))))
}

#[1] "number of individuals with greater than 10% missing data.....  33"
#[1] "number of individuals with greater than 20% missing data.....  1"
#[1] "number of individuals with greater than 30% missing data.....  0"
#[1] "number of individuals with greater than 40% missing data.....  0"
#[1] "number of individuals with greater than 50% missing data.....  0"
#[1] "number of individuals with greater than 60% missing data.....  0"
#[1] "number of individuals with greater than 70% missing data.....  0"
#[1] "number of individuals with greater than 80% missing data.....  0"
#[1] "number of individuals with greater than 90% missing data.....  0"

vcftools --vcf chdo_ULTIMO.vcf --012
vcftools --vcf chdo_ULTIMO.vcf --depth --out depth_ultima


After filtering, kept 268 Individuals, 18685 Sites, 27.49 X

vcftools --gzvcf chdo_ULTIMO.vcf --site-quality --out quality
vcftools --gzvcf chdo_ULTIMO.vcf --freq2 --out  --max-alleles 2
vcftools --gzvcf chdo_ULTIMO.vcf --depth --out meandepthind
vcftools --gzvcf chdo_ULTIMO.vcf --site-mean-depth --out meandepsite
vcftools --gzvcf chdo_ULTIMO.vcf --missing-indv --out missing
vcftools --gzvcf chdo_ULTIMO.vcf --missing-site --out missingsite
vcftools --gzvcf chdo_ULTIMO.vcf --het --out het
vcftools --gzvcf chdo_ULTIMO.vcf --counts --out count
```


### Genotype likelihood estimation

Updog and EBG:

```{r eval=FALSE}
1. Get sequencing error with updog, and then run ebg (following Trevor Faske)

### The vcf may have the following order: SpeciesPloidy_Pop_ID: CD2_AT_6 
Rscript updog_Rscript.R chdo_ultimo.vcf 4 2
Rscript updog_Rscript.R chdo_ultimo.vcf 6 2

require(readr)
library(data.table)
library(updog)
require(ggplot2)
library(vcfR)

## tetraploidies
updog4_out <- readRDS('/updog4_out.RDS')
seq_error4 <- updog4_out$snpdf$seq
hist(seq_error4)

updog4_snp <- as.character(updog4_out$snpdf$snp)
print(length(updog4_snp))
print(length(pos_ID))
print(updog4_snp[1:5])
print(pos_ID[1:5])
seq_err4_df <- data.frame(error=round(seq_error4,5))
write.table(seq_err4_df,'seq_error4n.txt',row.names=FALSE,col.names=FALSE)

## hexaploidies
updog6_out <- readRDS('/updog6_out.RDS')
seq_error6 <- updog6_out$snpdf$seq
hist(seq_error6)

updog6_snp <- as.character(updog6_out$snpdf$snp)
print(length(updog6_snp))
print(length(pos_ID))
print(updog6_snp[1:5])
print(pos_ID[1:5])
seq_err6_df <- data.frame(error=round(seq_error6,5))
write.table(seq_err6_df,'seq_error6n.txt',row.names=FALSE,col.names=FALSE)

##### EBG input files
vcf <- read.vcfR("chdo_ultimo.vcf", verbose = FALSE)

#get positions
chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
pos_ID <- paste(chrom,pos,sep = ':')

#Make pop_ID
indv <- colnames(ad)

Sp <- rep(NA,times=length(indv))
Ploidy <- rep(NA,times=length(indv))
Pop <- rep(NA,times=length(indv))
ID <- rep(NA,times=length(indv))
All <- rep(NA,times=length(indv))
for (i in 1:length(indv)){
  SpP <- unlist(strsplit(as.character(indv[i]),"_"))[1]
  Sp[i] <- gsub('\\d','',SpP,perl=TRUE)
  Ploidy[i] <-  gsub('(\\D)','',SpP,perl=TRUE)
  Pop[i] <- unlist(strsplit(as.character(indv[i]),"_"))[2]
  ID[i] <- unlist(strsplit(as.character(indv[i]),"_"))[3]
  All[i] <- as.character(indv[i])}
Pop_ID <- data.frame(Sp=Sp,Ploidy=Ploidy,Pop=Pop,ID=ID,All=All,
                     SpPloidy=paste0(Sp,Ploidy))
write.csv(Pop_ID,'Pop_ID.csv',row.names=FALSE)

#get AD
# ref, alt allele 
ad <- extract.gt(vcf, element = 'AD')
## check out PL and pos_ID
print(length(pos_ID))
print(pos_ID[1:10])
ad[1:5,1:5]

# get total, ref, and alt allele
tot_ad <- apply(ad, c(1,2), function(df) sum(as.numeric(unlist(strsplit(as.character(df),',')))))
ref_ad <- apply(ad, c(1,2), function(df) as.numeric(unlist(strsplit(as.character(df),','))[1]))
alt_ad <- apply(ad, c(1,2), function(df) as.numeric(unlist(strsplit(as.character(df),','))[2]))

print(ad[1:5,1:5])
print(dim(tot_ad))
print(ad[1:5,1:5])
print(ref_ad[1:5,1:5])
print(tot_ad[1:5,1:5])
print(alt_ad[1:5,1:5])

fwrite(tot_ad,'tot_ad.txt',quote=FALSE)
fwrite(ref_ad,'ref_ad.txt',quote=FALSE)
fwrite(alt_ad,'alt_ad.txt',quote=FALSE)

#tetraploids 4n
ploidy_index4 <- which(Pop_ID$Ploidy == 4)
tot_ad_4 <- tot_ad[,ploidy_index4]
ref_ad_4 <- ref_ad[,ploidy_index4]
alt_ad_4 <- alt_ad[,ploidy_index4]
Pop_ID_4 <- Pop_ID[ploidy_index4,]

print(length(ploidy_index4))
print(dim(tot_ad_4))
print(dim(ref_ad_4))
print(dim(alt_ad_4))

write.csv(Pop_ID_4,'Pop_ID_4.csv',row.names=FALSE)
fwrite(t(tot_ad_4),'tot_ad_4.txt',sep='\t',row.names=FALSE,col.names=FALSE)
fwrite(t(ref_ad_4),'ref_ad_4.txt',sep='\t',row.names=FALSE,col.names=FALSE)
fwrite(t(alt_ad_4),'alt_ad_4.txt',sep='\t',row.names=FALSE,col.names=FALSE)


###hexaploids 6n
ploidy_index6 <- which(Pop_ID$Ploidy == 6)
tot_ad_6 <- tot_ad[,ploidy_index6]
ref_ad_6 <- ref_ad[,ploidy_index6]
alt_ad_6 <- alt_ad[,ploidy_index6]
Pop_ID_6 <- Pop_ID[ploidy_index6,]

print(length(ploidy_index6))
print(dim(tot_ad_6))
print(dim(ref_ad_6))
print(dim(alt_ad_6))

write.csv(Pop_ID_6,'Pop_ID_6.csv',row.names=FALSE)
fwrite(t(tot_ad_6),'tot_ad_6.txt',sep='\t',row.names=FALSE,col.names=FALSE)
fwrite(t(ref_ad_6),'ref_ad_6.txt',sep='\t',row.names=FALSE,col.names=FALSE)
fwrite(t(alt_ad_6),'alt_ad_6.txt',sep='\t',row.names=FALSE,col.names=FALSE)


### Run ebg: 
ebg diseq -n 274 -l 22037 -t tot_ad_4.txt -a alt_ad_4.txt -e seqerror4.txt -p 4
ebg diseq -n 274 -l 22037 -t tot_ad_6.txt -a alt_ad_6.txt -e seqerror6.txt -p 6
```

### entropy: 1) all the populations as diploids:


```{r eval=FALSE}
perl vcf2mpglV1.3TLP.pl chdo_ULTIMO.vcf
perl gl2genestV1.3.pl chdo_ULTIMO.mpgl mean

# entropy 
require(readr)
require(MASS)
require(LEA)
require(ggplot2)
require(ggsci)
require(patchwork)

indv<-read.table("out.012.indv",sep="\t")

Sp <- rep(NA,times=nrow(indv))
Ploidy <- rep(NA,times=nrow(indv))
Pop <- rep(NA,times=nrow(indv))
ID <- rep(NA,times=nrow(indv))
All <- rep(NA,times=nrow(indv))
for (i in 1:nrow(indv)){
  SpP <- unlist(strsplit(as.character(indv$V1[i]),"_"))[1]
  Sp[i] <- gsub('\\d','',SpP,perl=TRUE)
  Ploidy[i] <-  gsub('(\\D)','',SpP,perl=TRUE)
  Pop[i] <- unlist(strsplit(as.character(indv$V1[i]),"_"))[2]
  ID[i] <- unlist(strsplit(as.character(indv$V1[i]),"_"))[3]
  All[i] <- as.character(indv$V1[i])}
Pop_ID <- data.frame(Sp=Sp,Ploidy=Ploidy,Pop=Pop,ID=ID,All=All,
                     SpPloidy=paste0(Sp,Ploidy))

print(head(Pop_ID))
write.csv(Pop_ID,"Pop_ID.csv",row.names = FALSE)

########### PCA_entropy
##row = loci, col = indv 
require(readr)
require(MASS)
require(LEA)
require(ggplot2)

PCA_entropy <- function(g){
  colmean<-apply(g,2,mean,na.rm=TRUE)
  normalize<-matrix(nrow = nrow(g),ncol=ncol(g))
  af<-colmean/2
  for (m in 1:length(af)){
    nr<-g[ ,m]-colmean[m]
    dn<-sqrt(af[m]*(1-af[m]))
    normalize[ ,m]<-nr/dn}
  
  normalize[is.na(normalize)]<-0
  method1<-prcomp(normalize, scale. = FALSE,center = FALSE)
  pve <- summary(method1)$importance[2,]
  print(pve[1:5])
  pca_df<- method1$x[,1:27]
  return(pca_df)} 

require(readr)
require(MASS)
require(LEA)
require(ggplot2)
g <- read.table("pntest_mean_chdo.txt", header=F)
Pop_ID_Sum <- read.csv("Pop_ID.csv")
pca_df <- PCA_entropy(t(g))

#Header for entropy:
Pop_ID <- read.csv("Pop_ID.csv")
Sp_Pop <- paste("CD",Pop_ID$Pop,sep="_")
Pop_ID <- paste(Pop_ID$Pop,Pop_ID$ID,sep="_")
Header <- data.frame(Sp_Pop,Pop_ID)
write.table(t(Header),'entropy_header.txt',sep = " ", quote = FALSE,row.names =FALSE)

dim(g)
#[1] loci_individuals

Pop_ID_Sum <- read.csv("Pop_ID.csv")
pca_df <- PCA_entropy(t(g)) 

k2<-kmeans(pca_df[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k3<-kmeans(pca_df[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k4<-kmeans(pca_df[,1:5],4,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k5<-kmeans(pca_df[,1:5],5,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k6<-kmeans(pca_df[,1:5],6,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k7<-kmeans(pca_df[,1:5],7,iter.max=10,nstart=10,algorithm="Hartigan-Wong")


ldak2<-lda(x=pca_df[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=pca_df[,1:5],grouping=k3$cluster,CV=TRUE)
ldak4<-lda(x=pca_df[,1:5],grouping=k4$cluster,CV=TRUE)
ldak5<-lda(x=pca_df[,1:5],grouping=k5$cluster,CV=TRUE)
ldak6<-lda(x=pca_df[,1:5],grouping=k6$cluster,CV=TRUE)
ldak7<-lda(x=pca_df[,1:5],grouping=k7$cluster,CV=TRUE)


write.table(round(ldak2$posterior,5),file="ldak2.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="ldak3.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak4$posterior,5),file="ldak4.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak5$posterior,5),file="ldak5.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak6$posterior,5),file="ldak6.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak7$posterior,5),file="ldak7.txt",quote=F,row.names=F,col.names=F)

*******
# in ponderosa:
#cat entropy_header.txt chdo_entropy.mpgl > entropy.mpgl

module load entropy/1.2
entropy -i entropy.mpgl -o chdo_k2.hdf5 -l 60000 -b 10000 -t 10 -s 20 -e .01 -k 2 -q ldak2.txt -m 1 -w 0 &> k2stdout.txt &
entropy -i entropy.mpgl -o chdo_k3.hdf5 -l 60000 -b 10000 -t 10 -s 20 -e .01 -k 3 -q ldak3.txt -m 1 -w 0 &> k3stdout.txt &
entropy -i entropy.mpgl -o chdo_k4.hdf5 -l 60000 -b 10000 -t 10 -s 20 -e .01 -k 4 -q ldak4.txt -m 1 -w 0 &> k4stdout.txt 
entropy -i entropy.mpgl -o chdo_k5.hdf5 -l 60000 -b 10000 -t 10 -s 20 -e .01 -k 5 -q ldak5.txt -m 1 -w 0 &> k5stdout.txt 
entropy -i entropy.mpgl -o chdo_k6.hdf5 -l 60000 -b 10000 -t 10 -s 20 -e .01 -k 6 -q ldak6.txt -m 1 -w 0 &> k6stdout.txt 
entropy -i entropy.mpgl -o chdo_k7.hdf5 -l 60000 -b 10000 -t 10 -s 20 -e .01 -k 7 -q ldak7.txt -m 1 -w 0 &> k7stdout.txt 
  
###Get the DICs values for each K value:

estpost.entropy chdo_k2.hdf5 -s 3 -p deviance > DIC_k2.txt
estpost.entropy chdo_k3.hdf5 -s 3 -p deviance > DIC_k3.txt
estpost.entropy chdo_k4.hdf5 -s 3 -p deviance > DIC_k4.txt*
estpost.entropy chdo_k5.hdf5 -s 3 -p deviance > DIC_k5.txt*
estpost.entropy chdo_k6.hdf5 -s 3 -p deviance > DIC_k6.txt*
estpost.entropy chdo_k7.hdf5 -s 3 -p deviance > DIC_k7.txt*

###Generate entropy q files:

#estpost.entropy chdo_k2.hdf5 -p q -s 0 -o q2.txt
#estpost.entropy chdo_k3.hdf5 -p q -s 0 -o q3.txt
#estpost.entropy chdo_k4.hdf5 -p q -s 0 -o q4.txt
#estpost.entropy chdo_k5.hdf5 -p q -s 0 -o q5.txt
#estpost.entropy chdo_k6.hdf5 -p q -s 0 -o q6.txt
#estpost.entropy basa_k7.hdf5 -p q -s 0 -o q7.txt

###Generate gprobs file:

#estpost.entropy  chdo_k2.hdf5 -p gprob -s 0 -o gprob.txt
```


### entropy: 2) with diploid likelihood from bcftools and hexa/tetraploid from ebg


```{r eval=FALSE}
require(readr)
library(data.table)
require(MASS)
require(LEA)
require(ggplot2)
library(vcfR)

setwd("/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/")

PL_4 <- data.frame(fread('4n/diseq-PL.txt',header=F,sep='\t'))
PL_6 <- data.frame(fread('6n/diseq-PL.txt',header=F,sep='\t'))

pntest_mean_2 <- data.frame(fread('pntest_mean_chdo_ULTIMO.txt'))
mpgl_2 <- data.frame(fread('chdo_ULTIMO.mpgl',header=FALSE, fill=TRUE),row.names=1)

Pop_ID_OG <- read.csv('Pop_ID_all.csv')
Pop_ID_4 <- read.csv('Pop_ID_4.csv')
Pop_ID_6 <- read.csv('Pop_ID_6.csv')

# REMOVE LAST COLUMN of PL_4/6 for some reason
PL_4 <- PL_4[,-c(ncol(PL_4))]
PL_6 <- PL_6[,-c(ncol(PL_6))]

print(dim(Pop_ID_OG))
print(dim(Pop_ID_4))

print(dim(pntest_mean_2))
print(pntest_mean_2[1:5,1:5])
print(dim(mpgl_2))
print(mpgl_2[1:5,1:5])
print(dim(PL_4))
print(PL_4[1:5,(ncol(PL_4)-5):ncol(PL_4)])
print(dim(PL_6))
print(PL_6[1:5,(ncol(PL_6)-5):ncol(PL_6)])

mpgl_4 <- apply(PL_4, c(1,2), function(df) gsub(',',' ',df,fixed=TRUE))
mpgl_6 <- apply(PL_6, c(1,2), function(df) gsub(',',' ',df,fixed=TRUE))

#check pl and mggl
print(PL_4[1:5,1:5])
print(mpgl_4[1:5,1:5])
print(PL_6[1:5,1:5])
print(mpgl_6[1:5,1:5])

rownames(mpgl_4) <- rownames(mpgl_2)
print(dim(mpgl_4))
fwrite(mpgl_4,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/mpgl_4.txt',quote=F,col.names = F,sep = ' ') 

rownames(mpgl_6) <- rownames(mpgl_2)
print(dim(mpgl_6))
fwrite(mpgl_6,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/mpgl_6.txt',quote=F,col.names = F,sep = ' ')

### pntest_mean ploidy###
mean_pl <- function(GP){
  gps <- as.numeric(unlist(strsplit(as.character(GP),split = ',')))
  ploidy <- length(gps)-1
  if (sum(gps) == 0){
    mean <- 0
  }else{
    gps <- 10 ^ (gps/-10)
    mean <- sum(gps*(0:ploidy))/sum(gps)
    return(mean) 
  }
}

system.time(pntest_mean_4 <- apply(PL_4,1:2,mean_pl))
system.time(pntest_mean_6 <- apply(PL_6,1:2,mean_pl))

rownames(pntest_mean_4) <- rownames(mpgl_2)
print(dim(pntest_mean_4))
fwrite(pntest_mean_4,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/pntest_mean_4.txt',quote=F,sep = ' ')

rownames(pntest_mean_6) <- rownames(mpgl_2)

print(dim(pntest_mean_6))
fwrite(pntest_mean_6,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/pntest_mean_6.txt',quote=F,sep = ' ')

#make sure ploidy makes sense
print(apply(pntest_mean_4,2,function(df) range(df, na.rm=TRUE)))[1:5]
print(apply(pntest_mean_6,2,function(df) range(df, na.rm=TRUE)))[1:5]

print(dim(pntest_mean_4))
print(pntest_mean_4[1:5,1:4])

print(dim(pntest_mean_6))
print(pntest_mean_6[1:5,1:4])

Pop_ID_OG <- read.csv('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/Pop_ID_all.csv')
Pop_ID_4 <- read.csv('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/Pop_ID_4.csv')
Pop_ID_6 <- read.csv('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/Pop_ID_6.csv')

pntest_mean_2 <- data.frame(fread('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/pntest_mean_chdo_ULTIMO.txt'))
mpgl_2 <- data.frame(fread('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/chdo_ULTIMO.mpgl',header=FALSE, fill = TRUE),row.names=1)

pntest_mean_4 <- data.frame(fread('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/pntest_mean_4.txt'))
mpgl_4 <- data.frame(fread('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/mpgl_4.txt',header=FALSE))

pntest_mean_6 <- data.frame(fread('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/pntest_mean_6.txt'))
mpgl_6 <- data.frame(fread('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/mpgl_6.txt',header=FALSE))

print(dim(Pop_ID_OG))
print(dim(pntest_mean_2))
print(dim(mpgl_2))

print(dim(Pop_ID_4))
print(dim(pntest_mean_4))
print(dim(mpgl_4))

print(dim(Pop_ID_6))
print(dim(pntest_mean_6))
print(dim(mpgl_6))

#extract only diploids for mpgl 
#first make a header for mpgl_2 so can select the right indv
names3 <- rep(Pop_ID_OG$All,each=3)
names(mpgl_2) <- names3
only2_names <- Pop_ID_OG$All[which(Pop_ID_OG$Ploidy == 2)]
index_mpgl2 <- which(names(mpgl_2) %in% only2_names)
mpgl_2only <- mpgl_2[,index_mpgl2]

#only diploids for rest
pntest_mean_2only <- pntest_mean_2[,which(Pop_ID_OG$Ploidy == 2)]
Pop_ID_2only <- Pop_ID_OG[which(Pop_ID_OG$Ploidy == 2),]

print(dim(Pop_ID_2only))
print(dim(pntest_mean_2only))
print(dim(mpgl_2only))

pntest_mean_ebgAll <- cbind(pntest_mean_2only,pntest_mean_4,pntest_mean_6)

mpgl_ebgAll <- cbind(mpgl_2only,mpgl_4,mpgl_6)
Pop_ID_ebgAll <- rbind(Pop_ID_2only,Pop_ID_4,Pop_ID_6)

print(dim(Pop_ID_ebgAll))
print(dim(pntest_mean_ebgAll))
print(dim(mpgl_ebgAll))

#test right dim for mpgl
print((nrow(Pop_ID_2only)*3) + (nrow(Pop_ID_4)*5) + (nrow(Pop_ID_6)*7))
print(pntest_mean_ebgAll[1:5,1:5])
print(rownames(mpgl_2)[1:10])
rownames(mpgl_ebgAll) <- rownames(mpgl_2)

rownames(mpgl_ebgAll) <- rownames(mpgl_2)
fwrite(mpgl_ebgAll,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/mpgl_ebgAll.txt',quote=F,col.names = F,row.names=TRUE,sep = ' ') 
#change pntest NAs to 0 
fwrite(pntest_mean_ebgAll,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/pntest_mean_ebgAll.txt',na=0,col.names = F,sep = ' ') 
write.csv(Pop_ID_ebgAll,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/Pop_ID_ebgAll.csv',row.names=F)

## create ldak files: accounting for ploidy

Pop_ID <- read.csv("/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/Pop_ID_ebgAll.csv")
pntest_mean <- fread("/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/pntest_mean_ebgAll.txt",header=F, data.table=F)
g <- t(pntest_mean)

###make sure all NA to 0 
#g[is.na(g)]<-0
print(length(which(is.na(g))))
print(dim(g))
g[1:5,1:5]

g[which(g == 0)] <- NA
print(length(which(is.na(g))))

#### create pl_list ####

pl_list <- Pop_ID$Ploidy

#### get stats on g and ploidy ####
nind <- nrow(g)
nloci <- ncol(g)
n2_index <- which(pl_list == 2)
n4_index <- which(pl_list == 4)
n6_index <- which(pl_list == 6)

#### check distributions  / range ###
g2 <- g[n2_index,]
g4 <- g[n4_index,]
g6 <- g[n6_index,]
pl_list <- Pop_ID$Ploidy

#### get stats on g and ploidy ####

nind <- nrow(g)
nloci <- ncol(g)
n2_index <- which(pl_list == 2)
n4_index <- which(pl_list == 4)
n6_index <- which(pl_list == 6)

#### check distributions  / range ###
g2 <- g[n2_index,]
g4 <- g[n4_index,]
g6 <- g[n6_index,]

#### scale withing ####
g2_scale <- scale(g2)
g4_scale <- scale(g4)
g6_scale <- scale(g6)
g_z <- rbind(g2_scale,g4_scale,g6_scale)

### if NAs
g_z[is.na(g_z)] <- 0

z_pca <- prcomp(as.matrix(g_z),scale. = F, center = F)
pve <- summary(z_pca)$importance[2,1:5]
print(summary(z_pca)$importance[2:3,1:5])

pca_df <- cbind(rbind(Pop_ID[n2_index,],Pop_ID[n4_index,],Pop_ID[n6_index,]),
                z_pca$x[,1:5])

group.colors   <- c(AI = "#F6EDC3", AO ="#E9967A", BB ="#66FFFF", CC = "#66CC33",CL = "#8B0000",CP ="#FF3030",CT ="#FF9933",IP = "#8B3A62", LA = "#FF82AB", 
                    LS ="#FFC0CB",MC ="#FFE1FF", MM="#CDB5CD",MP ="#8B7B8B", OV ="#9B30FF",PM="#8470FF",RC ="#27408B",SA ="#87CEFA", 
                    SB="#104E8B",SF="#8DB6CD",SO="#6E7B8B",TR="#FFB90F",UD="#EE7600",WM="#FFF68F",ST="#8B7500")

ggplot(data = pca_df, aes(x=PC1,y=PC2,
                          fill=Pop,shape=as.character(Ploidy))) + 
  geom_point(colour='black',size = 4) +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + 
  ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
  scale_fill_manual(name='Population:',values = group.colors) +
  scale_shape_manual(name='Ploidy:',values = c(21,22,24)) + 
  guides(fill = guide_legend(override.aes=list(pch=21))) +
  theme_bw() + 
  theme(#legend.position = 'none',
    axis.text = element_text(size=11), 
    axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave('/Users/carolinaosunamascaro/Desktop/chdo_ultimo/ebg/pca_ploidy_chdo.pdf',height=8,width=5)

###
k2<-kmeans(z_pca$x[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k3<-kmeans(z_pca$x[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k4<-kmeans(z_pca$x[,1:5],4,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k5<-kmeans(z_pca$x[,1:5],5,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k6<-kmeans(z_pca$x[,1:5],6,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k7<-kmeans(z_pca$x[,1:5],7,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

ldak2<-lda(x=z_pca$x[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=z_pca$x[,1:5],grouping=k3$cluster,CV=TRUE)
ldak4<-lda(x=z_pca$x[,1:5],grouping=k4$cluster,CV=TRUE)
ldak5<-lda(x=z_pca$x[,1:5],grouping=k5$cluster,CV=TRUE)
ldak6<-lda(x=z_pca$x[,1:5],grouping=k6$cluster,CV=TRUE)
ldak7<-lda(x=z_pca$x[,1:5],grouping=k7$cluster,CV=TRUE)

write.table(round(ldak2$posterior,5),file="/Users/carolinaosunamascaro/Desktop/chdo_ultimo/entropy/ldak2.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="/Users/carolinaosunamascaro/Desktop/chdo_ultimo/entropy/ldak3.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak4$posterior,5),file="/Users/carolinaosunamascaro/Desktop/chdo_ultimo/entropy/ldak4.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak5$posterior,5),file="/Users/carolinaosunamascaro/Desktop/chdo_ultimo/entropy/ldak5.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak6$posterior,5),file="/Users/carolinaosunamascaro/Desktop/chdo_ultimo/entropy/ldak6.txt",quote=F,row.names=F,col.names=F)                                      
write.table(round(ldak7$posterior,5),file="/Users/carolinaosunamascaro/Desktop/chdo_ultimo/entropy/ldak7.txt",quote=F,row.names=F,col.names=F)

ploidy_inds <- data.frame(Ploidy=Pop_ID$Ploidy)
write.table(ploidy_inds,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/entropy/ploidy_inds.txt',quote=F,row.names=F,col.names=F)

#### entropy headers

######### create entropy header ####

Pop_ID_list <- Pop_ID$All
Header <- data.frame(dims = NA,Pop_ID_list)
dim(pntest_mean)
df <- t(Header)
dims <- paste(dim(pntest_mean)[2],dim(pntest_mean)[1],sep = " ")
df[1,1] <- dims

write.table(df,'/Users/carolinaosunamascaro/Desktop/chdo_ultimo/entropy/entropy_header.txt',sep = " ",na ="",
            quote = FALSE,row.names = FALSE,col.names = FALSE)

```
