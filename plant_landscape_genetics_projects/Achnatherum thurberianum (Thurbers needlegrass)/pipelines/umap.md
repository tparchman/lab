### UMAP (following J. Jahner and T. Faske)

```{r eval=FALSE}

library(data.table)
library(ggplot2)
library(ggsci)
library(umap)
library(LEA)
library(readr)
library(ggpubr)

## PCA following Trevor Faske
PCA_gen <- function(df_gen,num=10,tw=FALSE,tw_pvalue=0.01){
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
  print(pve[1:5])
  
### adjust number of PC axes ###
  if(nrow(df_gen) < num){
    num <- nrow(df_gen)}
  
#### Tracy Widom, PC axes ####
  if(tw == TRUE){
    cat('\nRunning Tracy Widom test....\n\n')
    write.lfmm(normalize, "temp.lfmm")
    pca_tw <- pca("temp.lfmm",center=FALSE)
    tw <- tracy.widom(pca_tw)
    tw_sign <- tw$pvalues[tw$pvalues<=tw_pvalue]
    cat("\nNumber of TW sig. PC axes: ",length(tw_sign),'\n\n')
    num = length(tw_sign)
    unlink('temp.lfmm') }
  pca_X<- method1$x[,1:num]
  pca_X <- as.data.frame(pca_X)
  pca_out <- list("pca_df"=pca_X,"pve"=pve)
  return(pca_out)}

##############
g <- fread('gprob.txt',sep=',',data.table = F)
Pop_ID_Sum <- read.csv("Pop_ID.csv")

##### Run PCA #### 
pca_out <- PCA_gen(g,tw=TRUE)
pve <- pca_out$pve[1:5]
pve
ncol(pca_out$pca_df)
pca_df <- pca_out$pca_df 
pca_df <- cbind(Pop_ID_Sum,pca_df)

## Plot
col17 <- pal_d3(palette='category20')(20)
PCA_plot_1vs2 <- ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=as.character(Pop))) +
geom_point(colour='black',size = 4,pch=21) + ggtitle("PCA") +
xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
scale_fill_manual(values = col17) + 
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
PCA_plot_1vs2

col17 <- pal_d3(palette='category20')(20)
PCA_plot_2vs3 <- ggplot(data = pca_df, aes(x=PC2,y=PC3,fill=as.character(Pop))) +
geom_point(colour='black',size = 4,pch=21) + ggtitle("PCA") +
xlab(paste("PC",2," (",pve[2]*100,"%)",sep="")) + ylab(paste("PC",3," (",pve[3]*100,"%)",sep=""))  +
scale_fill_manual(values = col17) + 
theme_bw() +
theme(legend.position = 'none', #removes legend
     plot.title = element_text(size = 18, colour="black"),
    axis.text = element_text(size=16),
   axis.title = element_text(size = 18, colour="black",face = "bold"),    panel.border = element_rect(size = 1.5, colour = "black"),
       legend.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
      legend.text = element_text(size=13),
      panel.grid.major = element_blank(),
     panel.grid.minor = element_blank())
 PCA_plot_2vs3

 col17 <- pal_d3(palette='category20')(20)
 PCA_plot_3vs4 <- ggplot(data = pca_df, aes(x=PC3,y=PC4,fill=as.character(Pop))) +
   geom_point(colour='black',size = 4,pch=21) + ggtitle("PCA") +
   xlab(paste("PC",3," (",pve[3]*100,"%)",sep="")) + ylab(paste("PC",4," (",pve[4]*100,"%)",sep=""))  +
   scale_fill_manual(values = col17) + 
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
 PCA_plot_3vs4
 
 all_plots <- ggarrange(PCA_plot_1vs2,PCA_plot_2vs3,PCA_plot_3vs4, align='hv')
 all_plots
 ggsave("PCA_entropy.pdf",all_plots,height=8,width = 12,units = 'in')

 #### UMAP
custom.settings = umap.defaults
custom.settings$min_dist = 0.8
custom.settings$n_neighbors = 14

umap_g <- as.data.frame(umap(g,config = custom.settings)$layout )
names(umap_g) <- c('layout1','layout2')
umap_g <- cbind(Pop_ID_Sum,umap_g)

#### UMAP gprob ####
col17 <- pal_d3(palette='category20')(18)
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
 
