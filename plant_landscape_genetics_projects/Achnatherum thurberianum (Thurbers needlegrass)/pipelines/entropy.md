
---
---
#### Entropy (following T.Faske and J.Jahner notes)

Preparing files:

1. Generate mpgl file:

Needs vcf2mpglV1.3TLP.pl
```{r eval=FALSE}
perl /working/mascaro/acth/entropy2/vcf2mpglV1.3TLP.pl variants_maf5_miss9_thin100_noBadInds.recode.vcf
```

2. Create individuals file:
```{r eval=FALSE}
vcftools --vcf variants_maf5_miss9_thin100_noBadInds.recode.vcf --missing-indv 
cut -f 1 out.imiss > acth_individuals.txt
sed "s/INDV/ind/" acth_individuals.txt | sed "s/aln_//g" | sed "s/.sorted.bam//g" > acth_goodheads.txt
```

3. Generate pntest file:

Needs gl2genestV1.3.pl 

```{r eval=FALSE}
perl /working/mascaro/acth/entropy2/gl2genestV1.3.pl variants_maf5_miss9_thin100_noBadInds.mpgl mean
mv variants_maf5_miss9_thin100_noBadInds.mpgl acth_entropy.mpgl
```

4. Generate populatios file:
```{r eval=FALSE}
cut -d "_" -f 1,2 acth_goodheads.txt > acth_pops.txt
```

5. Make genotype likelihood matrix:

```{r eval=FALSE}
R
read.table("pntest_variants_maf5_miss9_thin100_noBadInds.txt", header=F)->gl
read.table("acth_individuos.txt", header=T)->ids
read.table("acth_pops.txt", header=T)->pops
t(gl)->tgl
cbind(ids, pops, tgl)->tidsgl
write.table(tidsgl, file="acth_16048_2round.txt", sep=" ", row.names=F, col.names=F , quote=F)
```

6. PCA for entropy:
```{r eval=FALSE}
R
PCA_entropy <- function(g){
  colmean<-apply(g,2,mean,na.rm=TRUE)
  normalize<-matrix(nrow = nrow(g),ncol=ncol(g))
  af<-colmean/2
  for (m in 1:length(af)){
    nr<-g[ ,m]-colmean[m]
    dn<-sqrt(af[m]*(1-af[m]))
    normalize[ ,m]<-nr/dn
  }
 
  normalize[is.na(normalize)]<-0
  method1<-prcomp(normalize, scale. = FALSE,center = FALSE)
  pve <- summary(method1)$importance[2,]
  print(pve[1:5])
  pca_df<- method1$x[,1:27]
  return(pca_df)
} 

require(readr)
require(MASS)
require(LEA)
require(ggplot2)
g <- read.table("pntest_mean_variants_maf5_miss9_thin100_noBadInds.recode.txt", header=F)
Pop_ID_Sum <- read.csv("Pop_ID.csv")
pca_df <- PCA_entropy(t(g))


Header for entropy:
Pop_ID <- read.csv("Pop_ID.csv")
Sp_Pop <- paste("GM",Pop_ID$Pop,sep="_")
Pop_ID <- paste(Pop_ID$Pop,Pop_ID$ID,sep="_")
Header <- data.frame(Sp_Pop,Pop_ID)
write.table(t(Header),'entropy_header.txt',sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE)


7. Generate LDA files
library(MASS)

k2<-kmeans(pca_df[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k3<-kmeans(pca_df[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k4<-kmeans(pca_df[,1:5],4,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k5<-kmeans(pca_df[,1:5],5,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k6<-kmeans(pca_df[,1:5],6,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k7<-kmeans(pca_df[,1:5],7,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k8<-kmeans(pca_df[,1:5],8,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k9<-kmeans(pca_df[,1:5],9,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k10<-kmeans(pca_df[,1:5],10,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

ldak2<-lda(x=pca_df[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=pca_df[,1:5],grouping=k3$cluster,CV=TRUE)
ldak4<-lda(x=pca_df[,1:5],grouping=k4$cluster,CV=TRUE)
ldak5<-lda(x=pca_df[,1:5],grouping=k5$cluster,CV=TRUE)
ldak6<-lda(x=pca_df[,1:5],grouping=k6$cluster,CV=TRUE)
ldak7<-lda(x=pca_df[,1:5],grouping=k7$cluster,CV=TRUE)
ldak8<-lda(x=pca_df[,1:5],grouping=k8$cluster,CV=TRUE)
ldak9<-lda(x=pca_df[,1:5],grouping=k9$cluster,CV=TRUE)
ldak10<-lda(x=pca_df[,1:5],grouping=k10$cluster,CV=TRUE)

write.table(round(ldak2$posterior,5),file="ldak2.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="ldak3.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak4$posterior,5),file="ldak4.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak5$posterior,5),file="ldak5.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak6$posterior,5),file="ldak6.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak7$posterior,5),file="ldak7.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak8$posterior,5),file="ldak8.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak9$posterior,5),file="ldak9.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak10$posterior,5),file="ldak10.txt",quote=F,row.names=F,col.names=F)
```

8. Make mpgl input file for entropy:

```{r eval=FALSE}
grep "_" acth_ids_good_head.txt > acth_nohead.txt
perl create_entropy_top_2rows.pl  acth_nohead.txt
cat entropy_2rows.txt acth_entropy.mpgl > acth_entropy2.mpgl
```

9. Running entropy (K 2-10):

```{r eval=FALSE}

module load entropy/1.2

entropy -i acth_entropy2.mpgl -o acth_k2.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 2 -q ldak2.txt -m 1 -w 0 &> k2stdout.txt &

entropy -i acth_entropy2.mpgl -o acth_k3.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 3 -q ldak3.txt -m 1 -w 0 &> k3stdout.txt &

entropy -i acth_entropy2.mpgl -o acth_k4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 4 -q ldak4.txt -m 1 -w 0 &> k4stdout.txt &

entropy -i acth_entropy2.mpgl -o acth_k5.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 5 -q ldak5.txt -m 1 -w 0 &> k5stdout.txt &

entropy -i acth_entropy2.mpgl -o acth_k6.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 6 -q ldak6.txt -m 1 -w 0 &> k6stdout.txt &

entropy -i acth_entropy2.mpgl -o acth_k7.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak7.txt -m 1 -w 0 &> k7stdout.txt &

entropy -i acth_entropy2.mpgl -o acth_k8.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak8.txt -m 1 -w 0 &> k8stdout4.txt &

entropy -i acth_entropy2.mpgl -o acth_k9.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak9.txt -m 1 -w 0 &> k9stdout4.txt &

entropy -i acth_entropy2.mpgl -o acth_k10.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak10.txt -m 1 -w 0 &> k10stdout.txt &
```

Get the DICs values for each K value:
```{r eval=FALSE}
module load entropy/1.2

estpost.entropy acth_k2.hdf5 -s 3 -p deviance
estpost.entropy acth_k3.hdf5 -s 3 -p deviance
estpost.entropy acth_k4.hdf5 -s 3 -p deviance
estpost.entropy acth_k5.hdf5 -s 3 -p deviance
estpost.entropy acth_k6.hdf5 -s 3 -p deviance
estpost.entropy acth_k7.hdf5 -s 3 -p deviance
estpost.entropy acth_k8.hdf5 -s 3 -p deviance
estpost.entropy acth_k9.hdf5 -s 3 -p deviance
estpost.entropy acth_k10.hdf5 -s 3 -p deviance
```
Generate entropy q files:

```{r eval=FALSE}
module load entropy/1.2

estpost.entropy acth_k2.hdf5 -p q -s 0 -o q2.txt
estpost.entropy acth_k3.hdf5 -p q -s 0 -o q3.txt
estpost.entropy acth_k4.hdf5 -p q -s 0 -o q4.txt
estpost.entropy acth_k5.hdf5 -p q -s 0 -o q5.txt
estpost.entropy acth_k6.hdf5 -p q -s 0 -o q6.txt
estpost.entropy acth_k7.hdf5 -p q -s 0 -o q7.txt
estpost.entropy acth_k8.hdf5 -p q -s 0 -o q8.txt
estpost.entropy acth_k9.hdf5 -p q -s 0 -o q9.txt
estpost.entropy acth_k10.hdf5 -p q -s 0 -o q10.txt
```

```{r eval=FALSE}
Generate gprobs file:

module load entropy/1.2

estpost.entropy  acth_k4.hdf5 -p gprob -s 0 -o gprob.txt
```
