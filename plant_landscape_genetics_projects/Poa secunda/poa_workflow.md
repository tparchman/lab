
Variant calling:

```{r eval=FALSE}
denovo_map.pl --samples /archive/parchman_lab/rawdata_to_backup/GSAF_11_20_bc_parse/POSE -T 5 -O /working/mascaro/pose/pop_map -m 3 -M 2 -n 2 -S -b 1
populations -b 1 -P /working/mascaro/pose -t 12 -M /working/mascaro/pose/pop_map --max_obs_het 0.65 -r 0.80 --vcf

```

Variant filtering:

```{r eval=FALSE}

vcftools --remove-indels --max-missing 0.4 --missing-indv --min-alleles 2 --max-alleles 2 --maf 0.01 --thin 5 --remove-filtered-all --recode --recode-INFO-all --vcf batch_1.vcf --out Poa_filter3

mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv

vcftools --vcf Poa_filter3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out variants_clean.recode.vcf

vcftools --vcf variants_clean.recode.vcf.recode.vcf --out Poa_filtered_final --remove-filtered-all --maf 0.01 --max-missing 0.4 --recode --thin 5


### 
After filtering, kept 260 out of 260 Individuals
After filtering, kept 9453 out of a possible 10034 Sites
###
```

PCA following Trevor Faske:


```{r eval=FALSE}
library(data.table)
library(tidyverse)
library(ggsci)
library('wesanderson')
library(readr)
library(readxl)

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

###### MAKE PCA 

setwd("/home/caro/Escritorio/Poa secunda/Poa filtro/")
df012<-fread("poa_filtered_012.012",sep="\t", data.table=F) 
df012 <- df012[,-1]

Pop_ID_Sum <- read.csv('Pop_ID.csv')

############# PCA ######################

pca_out <- PCA_gen(df012,Pop_ID_Sum)
pve <- pca_out$pve[1:5]
pve
#    PC1     PC2     PC3     PC4     PC5 
#  0.03349 0.03157 0.02907 0.02179 0.01969
#
pca_df <- pca_out$pca_df[,1:5]
pca_df <- cbind(Pop_ID_Sum,pca_df)

write.csv(pca_df,'pca_df.csv',row.names = FALSE)

############## PLOT ##################

pca_df <- read.csv('pca_df.csv')
pve <- c(0.03349, 0.03157, 0.02907, 0.02179, 0.01969)

#### Pop ####

PCA1VS2 <- ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=as.character(Pop))) + 
  geom_point(pch=21,colour='black',size = 4)+ #ggtitle("PCA Poa secunda, 6540 loci") +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
  #scale_fill_npg(name='Species:') + 
  #scale_fill_d3(name='Pop:',palette = 'category20') + 
  theme_bw() + 
  theme(legend.position = 'none',
    axis.text = element_text(size=11), 
    axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
ggsave("/Poa secunda/Poa filtro/pca1vs2.pdf",PCA1VS2,height=8,width = 12,units = 'in')


PCA2VS3 <- ggplot(data = pca_df, aes(x=PC2,y=PC3,fill=as.character(Pop))) + 
  geom_point(pch=21,colour='black',size = 4)+ #ggtitle("PCA snail ddocent") +
  xlab(paste("PC",2," (",pve[2]*100,"%)",sep="")) + ylab(paste("PC",3," (",pve[3]*100,"%)",sep=""))  +
  #scale_fill_npg(name='Species:') + 
  #scale_fill_d3(name='Pop:',palette = 'category20') + 
  theme_bw() + 
  theme(legend.position = 'none',
    axis.text = element_text(size=11), 
    axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
ggsave("/Poa secunda/Poa filtro/pca2vs3.pdf",PCA2VS3,height=8,width = 12,units = 'in')


PCA3VS4 <- ggplot(data = pca_df, aes(x=PC3,y=PC4,fill=as.character(Pop))) + 
  geom_point(pch=21,colour='black',size = 4)+ #ggtitle("PCA snail ddocent") +
  xlab(paste("PC",3," (",pve[3]*100,"%)",sep="")) + ylab(paste("PC",4," (",pve[4]*100,"%)",sep=""))  +
  #scale_fill_npg(name='Species:') + 
  #scale_fill_d3(name='Pop:',palette = 'category20') + 
  theme_bw() + 
  theme(legend.position = 'none',
    axis.text = element_text(size=11), 
    axis.title = element_text(size = 13, colour="black",face = "bold",vjust = 1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
ggsave("/Poa secunda/Poa filtro/pca3vs4.pdf",PCA3VS4,height=8,width = 12,units = 'in')

all_plots <- ggarrange(PCA1VS2,PCA2VS3,PCA3VS4, align='hv')
ggsave("Poa secunda/pca todas/all.pdf",all_plots,height=8,width = 12,units = 'in')
```
