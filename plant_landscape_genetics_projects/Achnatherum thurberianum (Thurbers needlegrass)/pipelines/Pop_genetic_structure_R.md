
#### Population genetic analyses in R from a vcf file
```{r eval=FALSE}
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(igraph)
library(reshape2)
library(ggplot2) 
library(cowplot)
library(tidyr)

# Open and explore the vcf file
ACTH.VCF <- read.vcfR("file.vcf")
ACTH.VCF

# Add a file with the populations information (Sample and ID)
pop.data <- read.table("pop_map_.txt", sep = "\t", header = TRUE)
##Check that all the samples in the VCF and the population data frame are included
all(colnames(ACTH.VCF@gt)[-1] == pop.data$ID)
##Convert the vcfR object into a genlight object
gl.rubi <- vcfR2genlight(ACTH.VCF)
pop(gl.rubi) <- pop.data$Pop
ploidy(gl.rubi) <- 2
gl.rubi

# Create a pairwise genetic distance matrix for individuals or populations 
rubi.dist <- poppr::bitwise.dist(gl.rubi)

##Minimum spanning networks
rubi.msn <- poppr.msn(gl.rubi, rubi.dist, showplot = FALSE, include.ties = T)
node.size <- rep(2, times = nInd(gl.rubi))
names(node.size) <- indNames(gl.rubi)
#vertex.attributes(rubi.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(gl.rubi, rubi.msn , palette = brewer.pal(n = nPop(gl.rubi), name = "Dark2"), gadj = 7, gweight = 2)

# Distance tree
tree <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
## color the tips of the tree based on the population of origin of the samples
cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.rubi)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend(35,10,c("","",""),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("","",""), fill = cols, border = FALSE, bty = "n", cex = 2)
#axis(side = 1)
#title(xlab = "Genetic distance (proportion of loci that are different)")

# PCA
rubi.pca <- glPca(gl.rubi, nf = 3)
barplot(100*rubi.pca$eig/sum(rubi.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
rubi.pca.scores <- as.data.frame(rubi.pca$scores)
rubi.pca.scores$pop <- pop(gl.rubi)
set.seed(9)
p <- ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p

# DAPC
pnw.dapc <- dapc(gl.rubi, n.pca = 3, n.da = 2)
cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2")
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)
compoplot(pnw.dapc,col = cols, posi = 'top')
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.rubi)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
write.csv(dapc.results, "dapc")
head(dapc.results, n = 18)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

##K means clustering
ACTH.VCF
##Add a file with the populations information (Sample and ID)
pop.data <- read.table("pop_map_.txt", sep = "\t", header = TRUE)

##Check that all the samples in the VCF and the population data frame are included
all(colnames(ACTH.VCF@gt)[-1] == pop.data$ID)

##Convert the vcfR object into a genlight object
gl.rubi <- vcfR2genlight(ACTH.VCF)
pop(gl.rubi) <- pop.data$Pop
maxK <- 18
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl.rubi, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat}
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1

my_k <- 2:18
grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(gl.rubi, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl.rubi, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

###############
# Extract genotypes
genotypes<- extract.gt(ACTH.VCF,element = "GT",mask = FALSE,as.numeric = FALSE,
                       return.alleles = FALSE,IDtoRowNames = TRUE,extract = TRUE,convertNA = TRUE)
write.csv(genotypes, "genotypes.csv")

###############
library(vcfR)
library("adegenet")
library("hierfstat")
library("pegas")
library(adegenet)

ACTH.VCF <- read.vcfR("file.vcf")
ACTH.VCF
x <- vcfR2genlight(ACTH.VCF)
pops <- as.factor(c("populationnames"))
ploidy(x) <- 2

# Genetic differentiation
x.dist <- poppr::bitwise.dist(x)
myDiff <- genetic_diff(ACTH.VCF, pops, method = 'nei')
write.csv(myDiff, "geneticdiff.csv")

# Population structure (Fis, Fst)
my_genind <- vcfR2genind(ACTH.VCF)

# Population specific Fis:
myData.hfstat <- genind2hierfstat(my_genind, pop = pops)
basicstat <- basic.stats(myData.hfstat, diploid = TRUE, digits = 4) 
basicstat$Fis
write.csv(basicstat$Fis, "Fis.csv")

# Bootstrapping over loci of population's Fis
boot.ppfis(myData.hfstat)

# Population specific FSTs:
betas(myData.hfstat,nboot=100,lim=c(0.025,0.975),diploid=TRUE)

# Nei's Pairwise FSTs: 
x <- genet.dist(myData.hfstat,diploid=TRUE,method="Ds")# Nei’s standard genetic distance

# Bootstrapping over loci of pairwise Fst
boot.ppfst(myData.hfstat)
```
