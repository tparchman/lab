
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
```

4. Generate populatios file:
```{r eval=FALSE}
cut -d "_" -f 1,2 acth_goodheads.txt > acth_pops.txt
```

```{r eval=FALSE}
5. Make genotype likelihood matrix:

R
read.table("pntest_variants_maf5_miss9_thin100_noBadInds.txt", header=F)->gl
read.table("acth_individuos.txtt", header=T)->ids
read.table("acth_pops.txt", header=T)->pops
t(gl)->tgl
cbind(ids, pops, tgl)->tidsgl
write.table(tidsgl, file="acth_16048_2round.txt", sep=" ", row.names=F, col.names=F , quote=F)
```

6. Generate LDA files
```{r eval=FALSE}
R
g <- read.table("acth_16048_2round.txt", header=F)
names <- read.table("acth_individuos.txt", header=T)
pops <- read.table("acth_pops.txt", header=T)
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
k8<-kmeans(pcg$x[,1:5],7,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

ldak2<-lda(x=pcg$x[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=pcg$x[,1:5],grouping=k3$cluster,CV=TRUE)
ldak4<-lda(x=pcg$x[,1:5],grouping=k4$cluster,CV=TRUE)
ldak5<-lda(x=pcg$x[,1:5],grouping=k5$cluster,CV=TRUE)
ldak6<-lda(x=pcg$x[,1:5],grouping=k6$cluster,CV=TRUE)
ldak7<-lda(x=pcg$x[,1:5],grouping=k7$cluster,CV=TRUE)
ldak8<-lda(x=pcg$x[,1:5],grouping=k7$cluster,CV=TRUE)

write.table(round(ldak2$posterior,5),file="ldak2.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="ldak3.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak4$posterior,5),file="ldak4.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak5$posterior,5),file="ldak5.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak6$posterior,5),file="ldak6.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak7$posterior,5),file="ldak7.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak7$posterior,5),file="ldak8.txt",quote=F,row.names=F,col.names=F)
```

7. Make mpgl input file for entropy:
```{r eval=FALSE}
grep "_" acth_goodheads.txt > acth_nohead.txt
perl /working/mascaro/acth/entropy2/create_entropy_top_2rows.pl acth_nohead.txt 
mkdir entropy2
mv ldak* entropy2/
mv entropy_2rows.txt entropy/
cd entropy/
cat entropy_2rows.txt variants_maf5_miss9_thin100_noBadInds.mpgl > acth_entropy.mpgl

*modify the header with a text editor (e.g.,nano) adding the numer of individuals, number of loci like that: 199 16048 1
```

8. Running entropy (K 2-10):

```{r eval=FALSE}
module load entropy/1.2

entropy -i acth_entropy.mpgl -o acth_k2.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 2 -q ldak2.txt -m 1 -w 0 &> k2stdout_rep4.txt &

entropy -i acth_entropy.mpgl -o acth_k3.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 3 -q ldak3.txt -m 1 -w 0 &> k3stdout_rep4.txt &

entropy -i acth_entropy.mpgl -o acth_k4.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 4 -q ldak4.txt -m 1 -w 0 &> k4stdout_rep4.txt &

entropy -i acth_entropy.mpgl -o acth_k5.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 5 -q ldak5.txt -m 1 -w 0 &> k5stdout_rep4.txt &

entropy -i acth_entropy.mpgl -o acth_k6.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 6 -q ldak6.txt -m 1 -w 0 &> k6stdout_rep4.txt &

entropy -i acth_entropy.mpgl -o acth_k7.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak7.txt -m 1 -w 0 &> k7stdout_rep4.txt &

entropy -i acth_entropy.mpgl -o acth_k8.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak7.txt -m 1 -w 0 &> k7stdout_rep4.txt &

entropy -i acth_entropy.mpgl -o acth_k9.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak7.txt -m 1 -w 0 &> k7stdout_rep4.txt &

entropy -i acth_entropy.mpgl -o acth_k10.hdf5 -l 70000 -b 30000 -t 10 -s 20 -e .01 -k 7 -q ldak7.txt -m 1 -w 0 &> k7stdout_rep4.txt &
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
