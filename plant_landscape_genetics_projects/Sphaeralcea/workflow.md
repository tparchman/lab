### SNP calling 

Sph data: 


S. ambigua (AM)
S. angustifolia (AN)
S. coccinea (CO)
S. fendleri (FE)
S. grossularifolia (GR)
S. hastulata (HA)
S. incana (IN)
S. munroana (MU)
S. parviflora (PA)
S. unknown (XX)


##### Generating the reference with dDocent:

```{r eval=FALSE}
mkdir refOpt (place here 4 individuals per pop)
./ReferenceOpt.sh 7 10 7 10 SE 5 0.95 0.99 0.005

#####
R
library(ggplot2)
data.table <- read.table("kopt.data", header = FALSE, col.names= c("k1","k2","Similarity", "Contigs"))
data.table$K1K2 <- paste(data.table$k1, data.table$k2, sep=",")
df=data.frame(data.table)
df$K1K2 <- as.factor(df$K1K2)
p <- ggplot(df, aes(x=Similarity, y=Contigs, group=K1K2)) + scale_x_continuous(breaks=seq(0.8,0.98,0.01)) + geom_line(aes(colour = K1K2))
p
ggsave("kvalues.pdf",p,height=8,width = 10,units = 'in')
#####

Run dDocent on this subset with the correct assembly parameters (skipping mapping and snp calling)

./RefMapOpt.sh 4 6 4 6 0.975 SE 10
./compress.sh (compress intermediate files)

Move the reads and reference files to the main directory
Run dDocent on the full data set, skipping trimming, assembly, snp calling
````

##### SNP calling:


```{r eval=FALSE}
Run ./bwa_sort.sh
INDS=($(for i in /working/mascaro/sph/*.F.fq.gz; do echo $(basename ${i%.F.fq.gz*}); done))

module load bwa/0.7.8
module load samtools/1.3

for IND in ${INDS[@]};
do
	# declare variables
	REF=/working/mascaro/sph/reference.fasta
	FORWARD=/working/mascaro/sph/${IND}.F.fq.gz
	OUTPUT=/working/mascaro/sph/assembly/${IND}_sort.bam
	# then align and sort
	echo "Aligning $IND with bwa"
	bwa mem -M -t 10 $REF $FORWARD \
	 | samtools view -b | \
	samtools sort -T ${IND} > $OUTPUT

done

###bcftools
Run ./bcftools.sh
REF=/working/mascaro/sph/reference.fasta

module load bcftools/1.9
bcftools mpileup -a AD,DP,SP -Ou -f $REF \
./assembly/*_sort.bam | bcftools call -f GQ,GP \
-mO z -o ./sph.vcf.gz
````

##### Filtering:


```{r eval=FALSE}
#### 

# Statistics with Vcftools:
vcftools --gzvcf sph.vcf.gz --site-quality --out /working/mascaro/sph/assembly/filter/quality
vcftools --gzvcf sph.vcf.gz --freq2 --out /working/mascaro/sph/assembly/filter --max-alleles 2
vcftools --gzvcf sph.vcf.gz --depth --out /working/mascaro/sph/assembly/filter/meandepthind
vcftools --gzvcf sph.vcf.gz --site-mean-depth --out /working/mascaro/sph/assembly/filter/meandepsite
vcftools --gzvcf sph.vcf.gz --missing-indv --out /working/mascaro/sph/assembly/filter/missing
vcftools --gzvcf sph.vcf.gz --missing-site --out /working/mascaro/sph/assembly/filter/missingsite

##### Examining statistics in R
R
library(ggplot2)
library(tidyverse)

var_qual <- read_delim("quality.lqual", delim = "\t",col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
ggsave("quality.pdf",a,height=8,width = 10,units = 'in')

var_depth <- read_delim("meandepthind.idepth", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip =1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_depth$mean_depth)
a + theme_light() + xlim(0, 100)
ggsave("meandepth_ind.pdf",a,height=8,width = 10,units = 'in')

var_miss <- read_delim("missingsite.lmiss", delim = "\t",col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss","fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)
ggsave("missing.pdf",a,height=8,width = 10,units = 'in')

ind_depth <- read_delim("meandepsite.ldepth.mean", delim = "\t", col_names = c("ind", "nsites", "depth"), skip =1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
ggsave("meandepth_site.pdf",a,height=8,width = 10,units = 'in')

ind_miss  <- read_delim("missing.imiss", delim = "\t",col_names = c("ind", "ndata", "nfiltered", "nmiss","fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
ggsave("missing.pdf",a,height=8,width = 10,units = 'in')
#######


# Filtering: 
Parameters as interpreted:
	--vcf sph_clean.recode.vcf
	--maf 0.05
	--thin 100
	--max-missing 0.8
	--out sph_maf005_thin100_maxmissing08
	--recode
	--remove-filtered-all



# Statistics with Vcftools:
vcftools --gzvcf sph_filtered_final.recode.vcf --site-quality --out quality
vcftools --gzvcf sph_filtered_final.recode.vcf --freq2 --out $OUT --max-alleles 2
vcftools --gzvcf sph_filtered_final.recode.vcf --depth --out meandepthind
vcftools --gzvcf sph_filtered_final.recode.vcf --site-mean-depth --out meandepsite
vcftools --gzvcf sph_filtered_final.recode.vcf --missing-indv --out missing
vcftools --gzvcf sph_filtered_final.recode.vcf --missing-site --out missingsite

After filtering, kept 560 out of 560 Individuals
After filtering, kept 19300 out of a possible 82583 Sites
```


##### Genetic diversity 

Genetic diversity statistics:

```{r eval=FALSE}
library(vcfR)
library("adegenet")
library("hierfstat")
library("pegas")
sph.vcf <- read.vcfR("sph_final.vcf")
sph.vcf

##Fis, Fst
my_genind <- vcfR2genind(sph.vcf)
x<- my_genind 
pops <- as.factor(c(pops))

#Population specific Fis:
myData.hfstat <- genind2hierfstat(my_genind, pop = pops)
basicstat <- basic.stats(myData.hfstat, diploid = TRUE, digits = 4) 
basicstat$Fis
write.csv(basicstat$Fis, "Fis.csv")

##Bootstrapping over loci of population's Fis
boot.ppfis(myData.hfstat)
#Nei's Pairwise FSTs: 
x <- genet.dist(myData.hfstat,diploid=TRUE,method="Ds")# Nei’s standard genetic distance
fst <- as.matrix(x)
write.table(fst, "Fst.txt")
##Bootstrapping over loci of pairwise Fst
#boot.ppfst(myData.hfstat)

basicstat$Ho #observed
write.csv(basicstat$Ho, "HO.csv")
basicstat$Hs #expected
write.csv(basicstat$Hs, "Hs.csv")
basicstat

###########################################
##vcftools Pi and TajimaD
#!/bin/sh
# .vcf file
# .pop file (unique names of pops, one per line)
# .map file (individual to population mapping file — 2 columns)

#cat sph.pop | while read line;
#do
#grep "$line" sph.map > $line.pop
#done

#for p in *.pop
#do
#vcftools --vcf sph.vcf --keep $p --site-pi --out $p
#done

#AMOVA:
library(vcfR)
library("adegenet")
library("hierfstat")
library("pegas")
library("poppr")
library("magrittr")
library(ape)

# vcf to genind
sph <- read.vcfR("sph.vcf")
my_genind <- vcfR2genind(sph)

### Give a data set with stratification (individual, populations, subpopulations..)
file.hier = read.table("ind_pop.txt", header=T)
strata(my_genind) = file.hier
my_genind

## amova with stratifications
amova.res.95 = poppr.amova(my_genind, ~pop/subpop, cutoff=0.95)
amova.res.95
write.table(amova.res.95$results, sep = ",", file = "results_amova.csv")
write.table(amova.res.95$componentsofcovariance, sep = ",", file = "Components_covariance.csv")
write.table(amova.res.95$statphi, sep = ",", file = "phi.csv")

## To test if populations are significantly different
amova.test.95 = randtest(amova.res.95)
amova.test.95
```

##### PCA, UMAP

```{r eval=FALSE}
#devtools::install_github("dkahle/ggmap")
library(ggmap)
library(devtools)
require(ggplot2)
require(ggsci)
require(ggrepel)
library(data.table)
library(patchwork)
options(ggrepel.max.overlaps = Inf)

setwd("/Proyectos/Postdoc/Sph/files/")
basacoord <- read.csv('coordinates2.csv', header=T, sep = ",")
names(basacoord) <- c('Pop','Lat','Long')

map <- get_stamenmap(bbox = c(left = -119.0, bottom = 30.00, right = -103.0, top = 49.00),
zoom=8,scale = 3,maptype = 'terrain-background')
##zoom controls resolution, 10 takes forever

ggmap(map) + geom_point(data = basacoord, aes(x = Long, y = Lat),size=3,pch=21,fill="white",col="black") +
xlab("Longitude") + ylab("Latitude") +
coord_map(xlim=c(-119.0,-103.0),ylim=c(30.00,49.00)) #selects the range for x and y

map <- ggmap(map) +
geom_point(data = basacoord, aes(x = Long, y = Lat),size=2,col="black",fill='white',pch=21) +
geom_label_repel(data = basacoord, aes(x = Long, y = Lat,label = Pop,fill=Pop, max.overlaps = Inf),
colour = "black", fontface = "bold") + scale_fill_manual(values=group.colors4) +
coord_map(xlim=c(-119.0,-103.0),ylim=c(30.00,49.00)) +
xlab("Longitude") + ylab("Latitude") + theme_bw() +
theme(legend.position = 'none',
axis.text = element_text(size=13),
axis.title = element_text(size = 15, colour="black",face = "bold", vjust = 1),
panel.border = element_rect(size = 1.5, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank())
map
ggsave("/map_pop_sph.pdf",height=8,width=5)
````


##### PCA and UMAP:

```{r eval=FALSE}
library(data.table)
library(ggplot2)
library(ggsci)
library(umap)
library(LEA)
library(readr)
library(ggpubr)

PCA_gen <- function(df_gen,indv,num=5){ #add ggplot, add tw, add # of output
  #pkgTest("ggplot2")
  df_gen <- apply(df_gen, 2, function(df) gsub(-1,NA,df,fixed=TRUE))
  df_gen <- apply(df_gen, 2, function(df) as.numeric(df))
  colmean<-apply(df_gen,2,mean,na.rm=TRUE)
  normalize<-matrix(nrow = nrow(df_gen),ncol=ncol(df_gen))
  af<-colmean/2
  for (m in 1:length(af)){
    nr<-df_gen[ ,m]-colmean[m]
    dn<-sqrt(af[m]*(1-af[m]))
    normalize[ ,m]<-nr/dn}
  
  normalize[is.na(normalize)]<-0
  method1<-prcomp(normalize, scale. = FALSE,center = FALSE)
  pve <- summary(method1)$importance[2,]
  print(pve[1:50])
  if(nrow(df_gen) < num){
    num <- nrow(df_gen)}
  
  pca_X<- method1$x[,1:num]
  pca_X <- as.data.frame(pca_X)
  pca_X$Pop <- indv$Pop
  pca_X$ID <- indv$ID
  pca_out <- list("pca_df"=pca_X,"pve"=pve)
  #print(PCA_fig(pca_out))
  return(pca_out)}

PCA_fig <- function(pca_out,fill="Pop"){ #add PCA setting or normal
  xlab <- paste("PCA1 (",round(pca_out$pve[1]*100,2),"%)",sep="")
  ylab <- paste("PCA2 (",round(pca_out$pve[2]*100,2),"%)",sep="")
  pca_df <- pca_out$pca_df
  filler <- as.character(pca_df[[fill]])
  ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=filler)) + 
    geom_point(pch=21,colour='black',size = 3)+ 
    xlab(xlab) + ylab(ylab) +
    theme_bw() + 
    theme(#legend.position = 'none',
      axis.text = element_text(size=13), 
      axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
      panel.border = element_rect(size = 1.5, colour = "black"),
      legend.text = element_text(size=13),
      legend.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())}

df012<-fread("out.012",sep="\t", data.table=F) 
df012 <- df012[,-1]
Pop_ID_Sum <- read.csv('Pop_ID.csv', sep = ";")


pca_out <- PCA_gen(df012,Pop_ID_Sum)
pve <- pca_out$pve[1:5]
pve

pca_df <- pca_out$pca_df[,1:5]
pca_df <- cbind(Pop_ID_Sum,pca_df)
write.csv(pca_df,'pca_df.csv',row.names = FALSE)

############## PLOT ##################
pca_df <- read.csv("pca_df.csv")
pve <- c(0.06978, 0.02266, 0.01435, 0.01289, 0.01046 ) 

PCA1VS2 <- ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=as.character(Pop), color=Pop)) + 
  geom_point(pch=21,colour='black',size = 5) +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))+
  scale_fill_manual(values=group.colors)+ theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=11), 
        axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("pca_spp.pdf",PCA1VS2,height=10,width = 10,units = 'in')

# PCA shapes per pop
PCA1VS2 <- ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=as.character(Pop), color=Pop)) + 
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))+
  geom_point(aes(shape=spp,size = 5)) +
  scale_color_manual(values=group.colors4)+
  scale_fill_manual(values=group.colors4)+ theme_bw() + 
  scale_shape_manual(values=c(15,7,9,18,17, 10,21,3,4,8,25)) +
  theme(legend.position = 'right',
        axis.text = element_text(size=11), 
        axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("pca_pop2.pdf",PCA1VS2,height=10,width = 10,units = 'in')

#### UMAP
custom.settings = umap.defaults
custom.settings$min_dist = 0.2
custom.settings$n_neighbors = 14

umap_g <- as.data.frame(umap(df012,config = custom.settings)$layout )
names(umap_g) <- c('layout1','layout2')
umap_g <- cbind(Pop_ID_Sum,umap_g)

#### UMAP gprob ####

umap_g_plot <- ggplot(data = umap_g, aes(x=layout1,y=layout2,fill=as.character(Pop))) +
  geom_point(colour='black',size = 5,pch=21) + #ggtitle("UMAP n_neighbors 14 min_dist 0.8") +
  xlab('Layout 1') + ylab('Layout 2') +
  scale_fill_manual(values=group.colors) + 
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size=11), 
        axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# UMAP shapes per pop
umap_g_plot <-ggplot(data = umap_g, aes(x=layout1,y=layout2,fill=as.character(Pop), color=Pop)) +
  #geom_point(colour='black',size = 5,pch=21) + #ggtitle("UMAP n_neighbors 14 min_dist 0.8") +
  xlab('Layout 1') + ylab('Layout 2') +
  geom_point(aes(shape=spp,size = 5)) +
  scale_color_manual(values=group.colors4)+
  scale_fill_manual(values=group.colors4)+ theme_bw() + 
  scale_shape_manual(values=c(15,7,9,18,17,10,21,3,4,8,25)) +
  theme(legend.position = 'right',
        axis.text = element_text(size=11), 
        axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
ggsave("umap_pop_label.pdf",umap_g_plot,height=10,width = 10,units = 'in')
```


##### Mean PCA, Procrustes:


```{r eval=FALSE}
setwd("/sph_pca_mean/")
library(vegan)
library(RColorBrewer)

mycols <- c("aliceblue", "darksalmon", "blanchedalmond", "lightcoral", "pelati", "hotpink3", "hotpink4", 
            "palevioletred", "pink","thistle1")
pcs_140 <- read.delim("pca_df.txt", header=FALSE, sep="")
	dim(pcs_140)

grouse <- read.delim("Summary_data.txt", header=TRUE)
	attach(grouse)
	dim(grouse)
	names(grouse)

######################################

quartz(width=4, height=12)
par(mar=c(4,4,2,1), mfrow=c(3,1))

## A) Total data fireworks

colorkey <- matrix(c("AM", mycols[1], "AN", mycols[2], "CO", mycols[3], "FE", mycols[4], "GR", mycols[5],"HA", mycols[6],"IN", mycols[7],"MU", mycols[8],"PA", mycols[9],"XX", mycols[10]), ncol=2, byrow=TRUE)

pcmeans <- data.frame(Lek=colorkey[,1], pc1=numeric(length(colorkey[,1])), pc2=numeric(length(colorkey[,1])))

for(i in 1:nrow(pcmeans)){
	pcmeans$pc1[i] <- mean(pcs_140[pcs_140[,1]==pcmeans$Lek[i], 2])
	pcmeans$pc2[i] <- mean(pcs_140[pcs_140[,1]==pcmeans$Lek[i], 3])	}

plot(-1*pcs_140[,3], -1*pcs_140[,2], type="n", cex.lab=1.25, cex.axis=1.25, las=1, xlab="PC 2 (3.376%)", ylab="PC 1 (1.843%)")
for(i in 1:nrow(pcmeans)){
	segments(-1*pcmeans$pc2[i], -1*pcmeans$pc1[i], -1*pcs_140[pcs_140[,1]==pcmeans$Lek[i], 3], -1*pcs_140[pcs_140[,1]==pcmeans$Lek[i], 2], col=colorkey[i,2], lwd=1.5)	}

points(-1*pcmeans$pc2[1], -1*pcmeans$pc1[1], pch=21, cex=2.5, lwd=1.25, bg=colorkey[1,2])
points(-1*pcmeans$pc2[2], -1*pcmeans$pc1[2], pch=21, cex=2.5, lwd=1.25, bg=colorkey[2,2])
points(-1*pcmeans$pc2[3], -1*pcmeans$pc1[3], pch=21, cex=2.5, lwd=1.25, bg=colorkey[3,2])
points(-1*pcmeans$pc2[4], -1*pcmeans$pc1[4], pch=21, cex=2.5, lwd=1.25, bg=colorkey[4,2])
points(-1*pcmeans$pc2[5], -1*pcmeans$pc1[5], pch=21, cex=2.5, lwd=1.25, bg=colorkey[5,2])
points(-1*pcmeans$pc2[6], -1*pcmeans$pc1[6], pch=21, cex=2.5, lwd=1.25, bg=colorkey[6,2])
points(-1*pcmeans$pc2[7], -1*pcmeans$pc1[7], pch=21, cex=2.5, lwd=1.25, bg=colorkey[7,2])
points(-1*pcmeans$pc2[8], -1*pcmeans$pc1[8], pch=22, cex=2.5, lwd=1.25, bg=colorkey[8,2])
points(-1*pcmeans$pc2[9], -1*pcmeans$pc1[9], pch=21, cex=2.5, lwd=1.25, bg=colorkey[9,2])
points(-1*pcmeans$pc2[10], -1*pcmeans$pc1[10], pch=21, cex=2.5, lwd=1.25, bg=colorkey[10,2])

mtext("A", adj=-0.15, cex=2, line=-1)

## B) Total data Procrustes

genes <- as.matrix(cbind(pcs_140[,2], pcs_140[,3]))
geography <- as.matrix(cbind(Long, Lat))

pro <- procrustes(geography, genes)
protest(geography, genes)

plotdata <- cbind(pro$Yrot, grouse$Lek)
mypalette <- brewer.pal(5, "Set1")

plot(pro, type="none", cex.axis=1.25, cex.lab=1.25, main="", xlab="Longitude", ylab="Latitude", las=1, ylim=c(-0.4, 0.3))
points(plotdata[plotdata[,3]=="AM", 1], plotdata[plotdata[,3]=="AM", 2], col=mypalette[1], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="AN", 1], plotdata[plotdata[,3]=="AN", 2], col=mypalette[2], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="CO", 1], plotdata[plotdata[,3]=="CO", 2], col=mypalette[3], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="FE", 1], plotdata[plotdata[,3]=="FE", 2], col=mypalette[4], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="GR", 1], plotdata[plotdata[,3]=="GR", 2], col=mypalette[5], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="HA", 1], plotdata[plotdata[,3]=="HA", 2], col=mypalette[6], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="IN", 1], plotdata[plotdata[,3]=="IN", 2], col=mypalette[7], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="MU", 1], plotdata[plotdata[,3]=="MU", 2], col=mypalette[8], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="PA", 1], plotdata[plotdata[,3]=="PA", 2], col=mypalette[9], pch=21, cex=1.25)
points(plotdata[plotdata[,3]=="XX", 1], plotdata[plotdata[,3]=="XX", 2], col=mypalette[10], pch=21, cex=1.25)

text(0.18, -0.35, labels="PC 1", cex=1.25)
text(0.4, 0.07, labels="PC 2", cex=1.25)
text(0.28, 0.3, labels="Procrustes Correlation = 0.695", cex=1.25)
mtext("B", adj=-0.15, cex=2, line=-1)
```


##### ADMIXTURE, piecharts


```{r eval=FALSE}
vcftools --vcf sph.vcf --plink-tped --out sph
plink --tped sph.tped --tfam sph.tfam --make-bed --out sph

for i in {2..10}
do
 ./admixture --cv sph.bed $i > log${i}.out
done

grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > binary.cv.error

# plotting
library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)

kdf10 <- read.csv("/q10.csv")
k10plot <-
  ggplot(kdf10, aes(factor(sampleID), prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(loc), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Populations", title = "K=10", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),axis.text.x = element_blank(),panel.grid = element_blank()) +
  scale_fill_gdocs(guide = FALSE)

k10plot
ggsave("k10.pdf",k8plot,height=8,width = 20,units = 'in')
```

##### Admixture piecharts:

```{r eval=FALSE}
## Admixture pie chart

BiocManager::install("LEA")
# Load packages
library(adegenet)
library(poppr)
library(LEA)
library(reshape2)
library(dplyr)
library(ggplot2)
library(rworldmap)
library(rworldxtra)
library(ggsn)
library(sf)
library(raster)
library(rgeos)
library(maps)
library(maptools)
library(grid)
library(miscTools)
library(stringr)
library(ggpubr)

qmatrix<-read.csv("/q10.csv",header=TRUE, sep = ";")
qmatrix
# Label column names of qmatrix
ncol(qmatrix)
head(qmatrix)

qlong = melt(qmatrix, id.vars=c("ID","Pop"))
head(qlong)

pal = colorRampPalette(c("red","blue", "green", "light blue"))
cols = pal(length(unique(qlong$variable)))

admix.bar = ggplot(data=qlong, aes(x=ID, y=value, fill=variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(~Pop, scales = "free", ncol = 2)+
  scale_fill_manual(values = cols)+
  ylab("Admixture proportion")+
  # xlab("Individual")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour="black", size=12),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
admix.bar
ggsave("Desktop/1.admixture_barplot.png", width=6, height=10, dpi=300)


# Calculate mean admixture proportions for each site
head(qmatrix)
clusters = grep("K", names(qmatrix)) # indexes of cluster columns
avg_admix = aggregate(qmatrix[, clusters], list(qmatrix$Pop), mean)

# Order alphabetically by site
avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix

# Convert dataframe from wide to long format
avg_admix = melt(avg_admix, id.vars = "Group.1")
head(avg_admix)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = cols)+
    theme_void()}

# Test function on one site
pie_charts(avg_admix, site = "AN", cols = cols)

# Subset data to reduce computation time
subsites = sort(c("AM", "AN", "CO", "FE", "GR", "HA", "IN", "MU", "PA", "XX"))

# Apply function to all sites using for loop
pies = list()
for (i in subsites){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = cols)}

# Prepare basemap
# Import csv file containing coordinates
coords = read.csv("coordinates.csv", sep = ";")

# Subset coordinates
coords = coords[coords$Pop %in% subsites, ]

# Order alphabetically by site
coords = coords[order(coords$Pop), ] 
coords

# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(coords$Pop)

# Set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-120.5, -109, 37, 49)
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Plot basemap
basemap = ggplot()+geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="grey",colour="black", size=0.5)+coord_quickmap(expand=F)+
  #lggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  #lggsn::scalebar(data = map.outline, dist = 50, dist_unit = "km", height = 0.01,
  #ltransform = TRUE, model = "WGS84" )+
                 #location = "bottomleft", anchor = c(x = -114, y = 38),
                 #st.bottom = FALSE, st.size = 4, st.dist = 0.015)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
    axis.text = element_text(colour="black", size=12),
    axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="lightsteelblue2"),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    panel.grid.minor = element_line(colour="grey90", size=0.5),
    panel.grid.major = element_line(colour="grey90", size=0.5),
    legend.text = element_text(size=12),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "top")
basemap

coord.list = list()
for (i in subsites){
  coord.list[[i]] = c(subset(coords, Pop == i)$Lon, subset(coords, Pop == i)$Lat)}
coord.list

# Define pie chart sizes
radius = 0.45

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(subsites)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)}

# Add layers to basemap
pie.map = basemap + pies.ac
pie.map
ggsave("Desktop/pie_charts_map.png", width = 8, height = 10, dpi = 300)

# Combine ggplots
ggarrange(admix.bar + theme(legend.position = "right") + labs(title = "Individual admixture proportions", tag = "A"),
          pie.map + labs(title = "Mean admixture proportions for each site", tag = "B"))
ggsave("Desktop/Admixture_bar_map.png", width = 30, height = 10, dpi = 300)
```

#####  IBD


```{r eval=FALSE}
library(ape)
library(ade4)
library(plyr)
library(vegan)
library(ggplot2)
library(dartR)

####IBD test
gen <- read.csv("pairwise.csv", header = TRUE)
Dgen <-as.matrix(gen)[, -1]
Dgen

coor <- read.csv("lat_lon.csv", header = TRUE)
latlong <- dist(cbind(coor$Lat, coor$Long))
Dgeo  <- as.matrix(latlong)[1:32, 1:32]
Dgeo

Dgen_2<-as.matrix(as.dist(Dgen))
gl.ibd(x = NULL,Dgen = Dgen_2,
       Dgeo = Dgeo,
       permutations = 999999,
       plot = TRUE)

######### Mantel statistic based on Spearman's rank correlation rho 

Mantel statistic based on Spearman's rank correlation rho 

Call:
mantel(xdis = Dgeo, ydis = Dgen, method = "spearman", permutations = 9999,      na.rm = TRUE) 

Mantel statistic r: 0.3322 
      Significance: 1e-04 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0702 0.0938 0.1145 0.1383 
Permutation: free
Number of permutations: 9999
```

#####  Phylogeny 

```{r eval=FALSE}
### vcf to phy
python vcf2phylip -i sph.vcf -b
### iqtree roostrap
./iqtree2 -s sph_filtered.phy -st DNA --model-joint 12.12 -b 1000 -t sph_filtered.phy.treefile -nt 15 --prefix sph_tree
```

```{r eval=FALSE}
### PyRAD 
 Parameters:
  
------ ipyrad params file (v.0.9.87)-------------------------------------------
test                           ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
/working/mascaro/sph/strong_filtering ## [1] [project_dir]: Project dir (made in curdir if not present)
                               ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                               ## [3] [barcodes_path]: Location of barcodes file
/working/parchman/SPH/rawfastqs/*fastq.gz ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference)
                               ## [6] [reference_sequence]: Location of reference sequence file
ddrad                          ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
TGCAG,                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                          ## [13] [maxdepth]: Max cluster depth within samples
0.9                            ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
0.05                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
0.05                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
20                             ## [21] [min_samples_locus]: Min # samples per locus for output
0.2                            ## [22] [max_SNPs_locus]: Max # SNPs per locus
8                              ## [23] [max_Indels_locus]: Max # of indels per locus
0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
p, s, l                        ## [27] [output_formats]: Output formats (see docs)
                               ## [28] [pop_assign_file]: Path to population assignment file
                               ## [29] [reference_as_filter]: Reads mapped to this reference are removed in step 3

##I've run a phylogeny with a subset of individuals from the phylip file obtained, and it looks like the former one-no ##species/pop clustering
##Now I'm running one with all the individuals (IQtree MFP+ASC) with the same phylip output of pyRAD
```

