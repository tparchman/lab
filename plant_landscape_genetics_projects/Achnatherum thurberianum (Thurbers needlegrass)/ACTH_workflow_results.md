---
title: "ACTH workflow and preliminary results"
author: ""
date: ""
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

######
  
## {.tabset}
### Variant calling and filtering

##### Running Stacks

```{r eval=FALSE}
module load stacks/1.46
module load gcc/4.8.5  

denovo_map.pl --samples /archive/parchman_lab/rawdata_to_backup/GSAF_11_20_bc_parse/ACTH/ -o /working/mascaro/ACTH/ -T 5 -O /home/working/ACTH/pop_map -m 3 -M 2 -n 2 -S -b 1
populations -b 1 -P /working/mascaro/ACTH -t 12 -M /working/mascaro/ACTH/pop_map --max_obs_het 0.65 -r 0.80 --vcf
```

##### Calculate allele frequency
```{r eval=FALSE}
vcftools --gzvcf VCF --freq2 --out $OUT --max-alleles 2
```
##### Calculate mean depth per individual
```{r eval=FALSE}
vcftools --gzvcf VCF --depth --out $OUT
```
##### Calculate mean depth per site
```{r eval=FALSE}
vcftools --gzvcf VCF --site-mean-depth --out $OUT
```
##### Calculate proportion of missing data per individual
```{r eval=FALSE}
vcftools --gzvcf VCF --missing-indv --out $OUT
```
##### Calculate proportion of missing data per site
```{r eval=FALSE}
vcftools --gzvcf VCF --missing-site --out $OUT
```

##### Examining statistics in R
```{r eval=FALSE}
library(tidyverse)

##Mean depth
var_depth <- read_delim("filtered.depth.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip =1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_depth$mean_depth)
a + theme_light() + xlim(0, 100)

##Variant missingness
var_miss <- read_delim("missing_site.txt.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss","fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)

##Minor alelle frecuencies
var_freq <- read_delim("frecuencies.txt.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_freq$maf)

##Mean depth per individual
ind_depth <- read_delim("depth.txt.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

##Proportion of missing data per individual
ind_miss  <- read_delim("missing_indiv.txt.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
```

##### vcf filtering with vcftools

```{r eval=FALSE}
vcftools --vcf batch_sin_.vcf --maf 0.03 --max-missing 0.8 --minQ 10 --min-meanDP 5 --max-meanDP 16 --minDP 5 --maxDP 16 --recode --out ACTH_filtered.vcf

mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv

vcftools --vcf ACTH_filtered.vcf --remove lowDP.indv --recode --recode-INFO-all --out variants_clean

vcftools --vcf variants_clean.recode.vcf --out ACTH_filtered_final --remove-filtered-all --maf 0.03 --max-missing 0.8 --recode --thin 5
```

##### __Final filtered vcf file__

```{r eval=FALSE}
5677 number of loci 
261 individuals
```

### Fis, Fst, Pi, TajimaD

Fis, Fst, Pi and TajimaD:

```{r eval=FALSE}
library(vcfR)
library("adegenet")
library("hierfstat")
library("pegas")
library(adegenet)
ACTH.VCF <- read.vcfR("ACTH.vcf")
ACTH.VCF
##Fis, Fst
my_genind <- vcfR2genind(ACTH.VCF)
x<- my_genind 
#Population specific Fis:
myData.hfstat <- genind2hierfstat(my_genind, pop = pops)
basicstat <- basic.stats(myData.hfstat, diploid = TRUE, digits = 4) 
basicstat$Fis
write.csv(basicstat$Fis, "Fis.csv")
##Bootstrapping over loci of population's Fis
boot.ppfis(myData.hfstat)
#Nei's Pairwise FSTs: 
x <- genet.dist(myData.hfstat,diploid=TRUE,method="Ds")# Nei’s standard genetic distance
##Bootstrapping over loci of pairwise Fst
boot.ppfst(myData.hfstat)

##vcftools Pi and TajimaD

vcftools --vcf ACTH.vcf --site-pi --positions pop.txt --out nucleotide_diversity
vcftools --vcf ACTH.vcf --out tajimasd --TajimaD 100
```


**Table 1.** Fis
```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='20%', fig.align="center"}
library(filenamer)
library(knitr)
library(kableExtra) 
library(dplyr)                                                         
data=read.csv(file = "/home/caro/Escritorio/fis.csv")
kable(data, caption = "") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE, italic = FALSE) %>% 
  #row_spec(7, underline = TRUE, monospace = TRUE) %>% 
  #row_spec(16, underline = TRUE, monospace = TRUE) %>%
  scroll_box(width = "1000px", height = "100px")
```



**Table 2.** Fst
```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='20%', fig.align="center"}
library(filenamer)
library(knitr)
library(kableExtra) 
library(dplyr)                                                         
data=read.csv(file = "/home/caro/Escritorio/fst.csv")
kable(data, caption = "") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE, italic = FALSE) %>% 
  #row_spec(7, underline = TRUE, monospace = TRUE) %>% 
  #row_spec(16, underline = TRUE, monospace = TRUE) %>%
  scroll_box(width = "1000px", height = "500px")
```


**Table 3.** Nucleotive diversity (Pi)
```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='20%', fig.align="center"}
library(filenamer)
library(knitr)
library(kableExtra) 
library(dplyr)                                                         
data=read.csv(file = "/home/caro/Escritorio/Values.csv")
kable(data, caption = "") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE, italic = FALSE) %>% 
  #row_spec(7, underline = TRUE, monospace = TRUE) %>% 
  #row_spec(16, underline = TRUE, monospace = TRUE) %>%
  scroll_box(width = "1000px", height = "100px")
```



**Table 4.** TajimaD
```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='20%', fig.align="center"}
library(filenamer)
library(knitr)
library(kableExtra) 
library(dplyr)                                                         
data=read.csv(file = "/home/caro/Escritorio/tajimaD.csv")
kable(data, caption = "") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE, italic = FALSE) %>% 
  #row_spec(7, underline = TRUE, monospace = TRUE) %>% 
  #row_spec(16, underline = TRUE, monospace = TRUE) %>%
  scroll_box(width = "1000px", height = "500px")
```


### PCA and UMAP

PCA following T. Faske:

```{r eval=FALSE}
library(data.table)
library(ggplot2)
library(ggsci)
library(umap)
library(LEA)
library(readr)
library(ggpubr)

##### PCA FUNCTION #####

##################### PCA_gen ##########################

#### PCA for 012 coded  vcf files
#### Following method in Patterson et al 2006
##### input files #####
### df_gen: genotypic data with individuals as rows and snps as coloumns.
###         Can include missing data. Either genotype probabilities or 012 format
### indv:  dataframe with rows corresponding to individuals in df_gen file 
###        Must have Pop and ID coloumn 
##### output files ####
### pca_out: 
### $ 'pca_df': dataframe with rows as individuals and coloumns as PC1-30, Pop, ID
### $ 'pve': list of proportion of variance explained for each PC

############# function ######################

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
    normalize[ ,m]<-nr/dn
  }
  
  normalize[is.na(normalize)]<-0
  method1<-prcomp(normalize, scale. = FALSE,center = FALSE)
  pve <- summary(method1)$importance[2,]
  print(pve[1:5])
  if(nrow(df_gen) < num){
    num <- nrow(df_gen)
  }
  
  pca_X<- method1$x[,1:num]
  pca_X <- as.data.frame(pca_X)
  pca_X$Pop <- indv$Pop
  pca_X$ID <- indv$ID
  pca_out <- list("pca_df"=pca_X,"pve"=pve)
  #print(PCA_fig(pca_out))
  return(pca_out)
}

######################## PCA_fig() ###################

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
      panel.grid.minor = element_blank())
}


######################################################################################
######################################################################################
#### START HERE #####
######################################################################################
######################################################################################
###### MAKE PCA 

df012<-fread("variantes012.012",sep="\t", data.table=F) 
df012 <- df012[,-1]
Pop_ID_Sum <- read.csv('Pop_ID.csv')

############# PCA ######################

pca_out <- PCA_gen(df012,Pop_ID_Sum)
pve <- pca_out$pve[1:5]
pve
#    PC1     PC2     PC3     PC4     PC5 
#0.15041 0.11894 0.11241 0.08590 0.05601 
pca_df <- pca_out$pca_df[,1:5]
pca_df <- cbind(Pop_ID_Sum,pca_df)

write.csv(pca_df,'pca_df.csv',row.names = FALSE)

############## PLOT ##################

pca_df <- read.csv("pca_df.csv")
pve <- c(0.15041, 0.11894, 0.11241, 0.08590, 0.05601)

#### Pop ####

PCA1VS2 <- ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=as.character(Pop))) + 
  geom_point(pch=21,colour='black',size = 4)+ #ggtitle("PCA ACTH") +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))+
  scale_fill_d3(name='Pop:',palette = 'category20') + 
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=11), 
        axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("pca1vs2_2.pdf",PCA1VS2,height=8,width = 12,units = 'in')


PCA2VS3 <- ggplot(data = pca_df, aes(x=PC2,y=PC3,fill=as.character(Pop))) + 
  geom_point(pch=21,colour='black',size = 4) +
  xlab(paste("PC",2," (",pve[2]*100,"%)",sep="")) + ylab(paste("PC",3," (",pve[3]*100,"%)",sep=""))+
  scale_fill_d3(name='Pop:',palette = 'category20') + 
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=11), 
        axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("pca2vs3.pdf",PCA2VS3,height=8,width = 12,units = 'in')


PCA3VS4 <- ggplot(data = pca_df, aes(x=PC3,y=PC4,fill=as.character(Pop))) + 
  geom_point(pch=21,colour='black',size = 4)+   #ggtitle("PCA") +
  xlab(paste("PC",3," (",pve[3]*100,"%)",sep="")) + ylab(paste("PC",4," (",pve[4]*100,"%)",sep=""))+
  scale_fill_d3(name='Pop:',palette = 'category20') + 
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=11), 
        axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("pca3vs4.pdf",PCA3VS4,height=8,width = 12,units = 'in')
```
  **-1.** PCA1 vs PCA2  
  **-2.** PCA2 vs PCA3  
  **-3.** PCA3 vs PCA4  
  
```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='100%', fig.align="center"}
library("EBImage")
setwd("/home/caro/Escritorio/figures_ACTH/PCA/")
net_name1 <- list.files(path="/home/caro/Escritorio/figures_ACTH/PCA/", pattern = "png")
display(readImage(net_name1, "PNG"), method="browser")
```

```{r eval=FALSE}
#### UMAP
custom.settings = umap.defaults
custom.settings$min_dist = 0.8
custom.settings$n_neighbors = 14

umap_g <- as.data.frame(umap(g,config = custom.settings)$layout )
names(umap_g) <- c('layout1','layout2')
umap_g <- cbind(Pop_ID_Sum,umap_g)

#### UMAP gprob ####
col17 <- pal_d3(palette='category20')(20)
umap_g_plot <- ggplot(data = umap_g, aes(x=layout1,y=layout2,fill=as.character(Pop))) +
  geom_point(colour='black',size = 4,pch=21) + ggtitle("UMAP n_neighbors 14 min_dist 0.8") +
  xlab('Layout 1') + ylab('Layout 2') +
  scale_fill_manual(values = col17) + 
  theme_bw() +
  theme(legend.position = 'none', 
        plot.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 18, colour="black",face = "bold"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        legend.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        legend.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
umap_g_plot
ggsave("umap_entropy.pdf",umap_g_plot,height=8,width = 8,units = 'in')
```

```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='100%', fig.align="center"}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage")
library("EBImage")
setwd("/home/caro/Escritorio/figures_ACTH/UMAP/")
net_name2 <- list.files(path="/home/caro/Escritorio/figures_ACTH/UMAP/", pattern = "png")
display(readImage(net_name2, "PNG"), method="browser")
```

### Isolation by distance

##### IBD:

```{r eval=FALSE}
#############
library(ape)
library(ade4)
library(plyr)
library(vegan)
library(ggplot2)
library(dartR)
####mantel test
gen <- read.csv("pairwise.csv")
Dgen <-as.matrix(gen)[, -1]
Dgen

coor <- read.csv("lat_lon.csv", header = TRUE)
latlong <- dist(cbind(coor$Lat, coor$Long))
Dgeo  <- as.matrix(latlong)[1:18, 1:18]
Dgeo

# mantel_vegan  <- mantel(Dgen, Dgeo, method = "spearman", permutations = 9999, na.rm = TRUE)
# mantel_vegan
# plot(Dgen, Dgeo)

Dgen_2<-as.matrix(as.dist(Dgen))
gl.ibd(x = NULL,Dgen = Dgen_2,
  Dgeo = Dgeo,
  permutations = 999999,
  plot = TRUE
)
```

__Results__:

```{r eval=FALSE}

#Mantel statistic based on Pearson's product-moment correlation 

Call:
vegan::mantel(xdis = Dgen, ydis = Dgeo, permutations = 999, na.rm = TRUE) 

Mantel statistic r: 0.4236 
      Significance: 0.001 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.137 0.177 0.226 0.266 
Permutation: free
Number of permutations: 999
```

```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='100%', fig.align="center"}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage")
library("EBImage")
setwd("/home/caro/Escritorio/figures_ACTH/IBD")
net_name3<- list.files(path="/home/caro/Escritorio/figures_ACTH/IBD", pattern = "png")
display(readImage(net_name3, "PNG"), method="browser")
```

### Admixture

Running admixture:

```{r eval=FALSE}
vcftools --vcf acth.vcf --plink-tped --out acth
plink --tped acth.tped --tfam acth.tfam --make-bed --out acth

for i in {2..15}
do
 ./admixture --cv acth.bed $i > log${i}.out
done

grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > binary.cv.error
Rscript admixture.R -p PH -i individual.populations.txt -k 8 -m 8 -l pop names
```

```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='100%', fig.align="center"}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage")
library("EBImage")
setwd("/home/caro/Escritorio/figures_ACTH/ADMIXTURE/")
net_name4<- list.files(path="/home/caro/Escritorio/figures_ACTH/ADMIXTURE/", pattern = "png")
display(readImage(net_name4, "PNG"), method="browser")
```

### Phylogenetic analyses

Convert vcf format to fasta using vcf2phylip.py:

```{r eval=FALSE}
python vcf2phylip.py --i acth.fasta --fasta
```

Trimming the alignment with trimAl:
  
```{r eval=FALSE}
./trimal -in acth.fasta -gappyout -out acth_fasta_trimmed.phy
```

Run iqtree:

```{r eval=FALSE}
./iqtree -s acth_fasta_trimmed.phy -nt 4 -st DNA -m MFP -bb 1000
```

```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='100%', fig.align="center"}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage")
library("EBImage")
setwd("/home/caro/Escritorio/figures_ACTH/Phylogenies1/")
net_name5<- list.files(path="/home/caro/Escritorio/figures_ACTH/Phylogenies1/", pattern = "png")
display(readImage(net_name5, "PNG"), method="browser")
```


```{r eval=FALSE}

library(phytools)
packageVersion("phytools")

# Latitude and longitude coordinates
coords<- read.csv("lat_long2.csv")
# Tree file
tree<-read.tree("acth_noname.fasta.tre")
rownames(coords)<-coords$Loc
coords<-coords[,-1]
library(mapdata)
library(viridis)
plot(tree)

# Drop tree (simplify the tree to get one sample per population, specify the tip to remove)
tree_drop <-drop.tip(tree, tip=c("AH_10","AH_11","AH_12","AH_13","AH_14","AH_2","AH_3","AH_4","AH_5","AH_6","AH_7","AH_8","AH_9","AS_10","AS_11","AS_12","AS_13","AS_14","AS_15","AS_2","AS_3","AS_4","AS_5","AS_6","AS_7","AS_8","AS_9","BM_10","BM_11","BM_12","BM_13","BM_14","BM_15","BM_16","BM_17","BM_18","BM_2","BM_3","BM_4","BM_5","BM_6","BM_7","BM_8","BM_9","BV_10","BV_11","BV_12","BV_13","BV_14","BV_15","BV_2","BV_3","BV_4","BV_5","BV_6","BV_7","BV_8","BV_9","DC_10","DC_11","DC_12","DC_13","DC_14","DC_15","DC_2","DC_3","DC_4","DC_5","DC_6","DC_7","DC_8","DC_9","DH_10","DH_11","DH_12","DH_13","DH_14","DH_2","DH_3","DH_4","DH_5","DH_6","DH_7","DH_8","DH_9","EW_10","EW_11","EW_12","EW_13","EW_14","EW_15","EW_2","EW_3","EW_4","EW_5","EW_6","EW_7","EW_8","EW_9","FR_10","FR_11","FR_12","FR_13","FR_14","FR_2","FR_3","FR_4","FR_5","FR_6","FR_7","FR_8","FR_9","GB_10","GB_11","GB_12","GB_13","GB_14","GB_2","GB_3","GB_4","GB_5","GB_6","GB_7","GB_8","GB_9","HO_10","HO_11","HO_12","HO_13","HO_14","HO_15","HO_2","HO_3","HO_4","HO_5","HO_6","HO_7","HO_8","HO_9","JC_10","JC_11","JC_12","JC_13","JC_14","JC_15","JC_2","JC_3","JC_4","JC_5","JC_6","JC_7","JC_8","JC_9","LV_10","LV_11","LV_12","LV_13","LV_14","LV_2","LV_3","LV_4","LV_5","LV_6","LV_7","LV_8","LV_9","MD_10","MD_11","MD_12","MD_13","MD_14","MD_15","MD_2","MD_3","MD_4","MD_5","MD_6","MD_7","MD_8","MD_9","PL_10","PL_11","PL_12","PL_13","PL_14","PL_15","PL_2","PL_3","PL_4","PL_5","PL_6","PL_7","PL_8","PL_9","PT_10","PT_11","PT_12","PT_2","PT_3","PT_4","PT_5","PT_6","PT_7","PT_8","PT_9","SC_10","SC_11","SC_12","SC_13","SC_14","SC_15","SC_2","SC_3","SC_4","SC_5","SC_6","SC_7","SC_8","SC_9","SS_10","SS_11","SS_12","SS_13","SS_14","SS_2","SS_3","SS_4","SS_5","SS_6","SS_7","SS_8","SS_9","VM_10","VM_11","VM_12","VM_2","VM_3","VM_4","VM_5","VM_6","VM_7", "VM_8", "VM_9"), trim.internal = TRUE, subtree = FALSE,
                 root.edge = 0,collapse.singles = TRUE,
                 interactive = FALSE)
tree_drop$tip.label
# Load the database "state" and specify which states would be pictured
obj<-phylo.to.map(tree_drop,coords,database="state",
                  regions=c("Nevada","California","Oregon"),plot=F)

plot(obj,direction="rightwards",colors=cols,xlim=c(-100,100),ylim=c(-100,100))
plot(obj,type="phylogram",asp=1)
```


```{r echo=FALSE, message=FALSE, warning=FALSE,out.width='100%', fig.align="center"}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage")
library("EBImage")
setwd("/home/caro/Escritorio/figures_ACTH/Phylogenies2/")
net_name6<- list.files(path="/home/caro/Escritorio/figures_ACTH/Phylogenies2/", pattern = "png")
display(readImage(net_name6, "PNG"), method="browser")
```

