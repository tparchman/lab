##########################################
## 1) load data and packages
##########################################

library(ape)
library(ecodist)
library(fossil)

ids <- read.delim("grouse2_ids_112.txt", header=TRUE)
	dim(ids)
	head(ids)

pops <- read.delim("grouse2_pops_112.txt", header=TRUE)
	dim(pops)
	head(pops)

gprobs <- read.csv("gprob2_rep0.txt", header=TRUE)
	dim(gprobs)
	gprobs[1:10,1:10]

gprobs_noname <- gprobs[,-c(1:2)]
	dim(gprobs_noname)
	gprobs_noname[1:10,1:10]


##################################################
## 5) umap connectivity
##################################################

library(umap)

pca <- prcomp(gprobs_noname, center=TRUE, scale=FALSE)
	summary(pca)

ninds <- dim(gprobs_noname)[1]
npcs <- 10
umap_pcs <- umap(pca$x[,1:npcs])

quartz(height=6, width=6)
par(mar=c(5,5,1,1))
plot(umap_pcs$layout[,1], umap_pcs$layout[,2], type="n", xlab="Layout 1", ylab="Layout 2", cex.lab=2, cex.axis=1.5, las=1)
for(i in 1:ninds)
	{
	for (j in 2:15)
		{
		line_weight <- 1.01 - (umap_pcs$knn$distances[i,j]/max(umap_pcs$knn$distances))
		segments(umap_pcs$layout[i,1], umap_pcs$layout[i,2], umap_pcs$layout[umap_pcs$knn$indexes[i,j],1], umap_pcs$layout[umap_pcs$knn$indexes[i,j],2], lwd=line_weight)
		}
	}
points(umap_pcs$layout[pops[,1]=="HA_HA", 1], umap_pcs$layout[pops[,1]=="HA_HA", 2], pch=21, bg="#c09700", cex=2)
points(umap_pcs$layout[pops[,1]=="MA_FA", 1], umap_pcs$layout[pops[,1]=="MA_FA", 2], pch=21, bg="#d10092", cex=2)
points(umap_pcs$layout[pops[,1]=="MA_NB", 1], umap_pcs$layout[pops[,1]=="MA_NB", 2], pch=22, bg="#d10092", cex=2)
points(umap_pcs$layout[pops[,1]=="MA_TL", 1], umap_pcs$layout[pops[,1]=="MA_TL", 2], pch=23, bg="#d10092", cex=2)
points(umap_pcs$layout[pops[,1]=="SH_HL", 1], umap_pcs$layout[pops[,1]=="SH_HL", 2], pch=21, bg="#0076da", cex=2)
points(umap_pcs$layout[pops[,1]=="SH_LS", 1], umap_pcs$layout[pops[,1]=="SH_LS", 2], pch=22, bg="#0076da", cex=2)
legend("topright", legend=c("HA_HA", "MA_FA", "MA_NB", "MA_TL", "SH_HL", "SH_LS"), pch=c(21,21,22,23,21,22), pt.bg=c("#c09700", "#d10092", "#d10092", "#d10092", "#0076da", "#0076da"), cex=1.5, pt.cex=2)



