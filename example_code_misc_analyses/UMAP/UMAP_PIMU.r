#### UMAP_PIMU.R
#### following Faske, Jahner, Diaz-Papkovich 2019 Plos Genetics
#### T. Parchman
#### 06/03/2021

library(data.table)
library(ggplot2)
library(ggsci)
library(umap)
library(LEA)
library(readr)
library(ggpubr)


################ FUNCTIONS ###############

##################### PCA_gen ##########################

#### PCA for 012 coded  vcf files
#### Following method in Patterson et al 2006

##### input files #####

### df_gen: genotypic data with individuals as rows and snps as coloumns.
###         Can include missing data. Either genotype probabilities or 012 format


##### output files ####

### pca_out: 

### $ 'pca_df': dataframe with rows as individuals and columns as PC1-X, Pop, ID
### $ 'pve': list of proportion of variance explained for each PC

############# function ######################

PCA_gen <- function(df_gen,num=10,tw=FALSE,tw_pvalue=0.01){
  
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
  
  ### adjust number of PC axes ###
  
  if(nrow(df_gen) < num){
    num <- nrow(df_gen)
  }
  
  #### Tracy Widom, PC axes ####
  if(tw == TRUE){
    cat('\nRunning Tracy Widom test....\n\n')
    write.lfmm(normalize, "temp.lfmm")
    pca_tw <- pca("temp.lfmm",center=FALSE)
    tw <- tracy.widom(pca_tw)
    tw_sign <- tw$pvalues[tw$pvalues<=tw_pvalue]
    cat("\nNumber of TW sig. PC axes: ",length(tw_sign),'\n\n')
    num = length(tw_sign)
    unlink('temp.lfmm')
  }
  
  pca_X<- method1$x[,1:num]
  
  pca_X <- as.data.frame(pca_X)
  
  pca_out <- list("pca_df"=pca_X,"pve"=pve)
  
  return(pca_out)
}

###### EXAMPLE: PIMU all ##########

#### setwd ####
setwd('~/Desktop/files/cal_serot_pines/calserpines_plotting_code/PIMU/UMAP/')

#### read in files ####
g <- fread('PM_gl_matrix_miss30_maf05_noBadInds_noHighCov_noParalogs_noWeird.recode.csv',sep=',',data.table = F)
g <- g[,-c(1:2)]

Pop_ID_Sum <- read.csv("PM_pop_ids.csv")

##### Run PCA #### 
pca_out <- PCA_gen(g,tw=TRUE)
pve <- pca_out$pve[1:5]
pve
#       PC1     PC2     PC3     PC4     PC5
#    0.07045 0.02406 0.01839 0.01187 0.01119

ncol(pca_out$pca_df) # 14, number of tw PC axes

pca_df <- pca_out$pca_df 
pca_df <- cbind(Pop_ID_Sum,pca_df)

#### UMAP ####

umap_g <- as.data.frame(umap(g)$layout)
names(umap_g) <- c('layout1','layout2')
umap_g <- cbind(Pop_ID_Sum,umap_g)

umap_tw_pcs <- as.data.frame(umap(pca_out$pca_df)$layout) #number of tw PC axes
names(umap_tw_pcs) <- c('layout1','layout2')
umap_tw_pcs <- cbind(Pop_ID_Sum,umap_tw_pcs)

umap_ten_pcs <- as.data.frame(umap(pca_out$pca_df[,1:10])$layout)
names(umap_ten_pcs) <- c('layout1','layout2')
umap_ten_pcs <- cbind(Pop_ID_Sum,umap_ten_pcs)


#### PC colour=color####

pve <- c(0.07045 0.02406 0.01839 0.01187 0.01119)

#col17 <- pal_d3(palette='category20')(20)[c(1:5,7,9:17,19,20)]
PCA_plot <- ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=as.character(Pop))) +
  geom_point(colour='black',size = 4,pch=21) + ggtitle("PCA") +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
# scale_fill_manual(values = col17) +
  scale_fill_d3(palette = 'category20') +
  theme_bw() +
theme(#legend.position = 'none', #removes legend
        plot.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 18, colour="black",face = "bold"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        legend.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        legend.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
PCA_plot

#### UMAP gprob ####

#col17 <- pal_d3(palette='category20')(20)[c(1:5,7,9:17,19,20)]
umap_g_plot <- ggplot(data = umap_g, aes(x=layout1,y=layout2,fill=as.character(Pop))) +
  geom_point(colour='black',size = 4,pch=21) + ggtitle("UMAP gprobs") +
  xlab('Layout 1') + ylab('Layout 2') +
#scale_fill_manual(values = col17) +
  scale_fill_d3(palette = 'category20') +
  theme_bw() +
  theme(legend.position = 'none', #removes legend
        plot.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 18, colour="black",face = "bold"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        legend.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        legend.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
umap_g_plot

#### UMAP tracy widom ####

#col17 <- pal_d3(palette='category20')(20)[c(1:5,7,9:17,19,20)]
umap_tw_plot <- ggplot(data = umap_tw_pcs, aes(x=layout1,y=layout2,fill=as.character(Pop))) +
  geom_point(colour='black',size = 4,pch=21) + ggtitle(paste0("UMAP tw, ",ncol(pca_out$pca_df)," PC axes")) +
  xlab('Layout 1') + ylab('Layout 2') +
# scale_fill_manual(values = col17) +
  scale_fill_d3(palette = 'category20') +
  theme_bw() +
  theme(legend.position = 'none', #removes legend
        plot.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 18, colour="black",face = "bold"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        legend.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        legend.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
umap_tw_plot

#### UMAP 10 pcs ####

#col17 <- pal_d3(palette='category20')(20)[c(1:5,7,9:17,19,20)]
umap_ten_plot <- ggplot(data = umap_ten_pcs, aes(x=layout1,y=layout2,fill=as.character(Pop))) +
  geom_point(colour='black',size = 4,pch=21) + ggtitle("UMAP 10 PC axes") +
  xlab('Layout 1') + ylab('Layout 2') +
# scale_fill_manual(values = col17) +
  scale_fill_d3(palette = 'category20') +
  theme_bw() +
  theme(legend.position = 'none', #removes legend
        plot.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 18, colour="black",face = "bold"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        legend.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        legend.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
umap_ten_plot

### Combine plots ####

all_plots <- ggarrange(PCA_plot,umap_tw_plot,umap_ten_plot,ncol=3)
all_plots

ggsave('PIMU_PCA_UMAP.pdf',all_plots,height=5,width = 20,units = 'in')
