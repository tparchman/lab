
### Projecting a phylogenetic tree onto a map:

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
tree_drop <-drop.tip(tree, tip=c("AH_10","AH_11","AH_12","AH_13","AH_14","AH_2","AH_3",	"AH_4",	"AH_5",	"AH_6",	"AH_7",	"AH_8",	"AH_9",	"AS_10","AS_11","AS_12","AS_13","AS_14","AS_15","AS_2",	"AS_3",	"AS_4",	"AS_5",	"AS_6",	"AS_7",	"AS_8",	"AS_9",	"BM_10",	"BM_11",	"BM_12",	"BM_13",	"BM_14",	"BM_15",	"BM_16",	"BM_17",	"BM_18",	"BM_2",	"BM_3",	"BM_4",	"BM_5",	"BM_6",	"BM_7",	"BM_8",	"BM_9",	"BV_10",	"BV_11",	"BV_12",	"BV_13",	"BV_14",	"BV_15",	"BV_2",	"BV_3",	"BV_4",	"BV_5",	"BV_6",	"BV_7",	"BV_8",	"BV_9",	"DC_10",	"DC_11",	"DC_12",	"DC_13",	"DC_14",	"DC_15",	"DC_2",	"DC_3",	"DC_4",	"DC_5",	"DC_6",	"DC_7",	"DC_8",	"DC_9",	"DH_10",	"DH_11",	"DH_12",	"DH_13",	"DH_14",	"DH_2",	"DH_3",	"DH_4",	"DH_5",	"DH_6",	"DH_7",	"DH_8",	"DH_9",	"EW_10",	"EW_11",	"EW_12",	"EW_13",	"EW_14",	"EW_15",	"EW_2",	"EW_3",	"EW_4",	"EW_5",	"EW_6",	"EW_7",	"EW_8",	"EW_9",	"FR_10",	"FR_11",	"FR_12",	"FR_13",	"FR_14",	"FR_2",	"FR_3",	"FR_4",	"FR_5",	"FR_6",	"FR_7",	"FR_8",	"FR_9",	"GB_10",	"GB_11",	"GB_12",	"GB_13",	"GB_14",	"GB_2",	"GB_3",	"GB_4",	"GB_5",	"GB_6",	"GB_7",	"GB_8",	"GB_9",	"HO_10",	"HO_11",	"HO_12",	"HO_13",	"HO_14",	"HO_15",	"HO_2",	"HO_3",	"HO_4",	"HO_5",	"HO_6",	"HO_7",	"HO_8",	"HO_9",	"JC_10",	"JC_11",	"JC_12",	"JC_13",	"JC_14",	"JC_15",	"JC_2",	"JC_3",	"JC_4",	"JC_5",	"JC_6",	"JC_7",	"JC_8",	"JC_9",	"LV_10",	"LV_11",	"LV_12",	"LV_13",	"LV_14",	"LV_2",	"LV_3",	"LV_4",	"LV_5",	"LV_6",	"LV_7",	"LV_8",	"LV_9",	"MD_10",	"MD_11",	"MD_12",	"MD_13",	"MD_14",	"MD_15",	"MD_2",	"MD_3",	"MD_4",	"MD_5",	"MD_6",	"MD_7",	"MD_8",	"MD_9",	"PL_10",	"PL_11",	"PL_12",	"PL_13",	"PL_14",	"PL_15",	"PL_2",	"PL_3",	"PL_4",	"PL_5",	"PL_6",	"PL_7",	"PL_8",	"PL_9",	"PT_10",	"PT_11",	"PT_12",	"PT_2",	"PT_3",	"PT_4",	"PT_5",	"PT_6",	"PT_7",	"PT_8",	"PT_9",	"SC_10",	"SC_11",	"SC_12",	"SC_13",	"SC_14",	"SC_15",	"SC_2",	"SC_3",	"SC_4",	"SC_5",	"SC_6",	"SC_7",	"SC_8",	"SC_9",	"SS_10",	"SS_11",	"SS_12",	"SS_13",	"SS_14",	"SS_2",	"SS_3",	"SS_4",	"SS_5",	"SS_6",	"SS_7",	"SS_8",	"SS_9",	"VM_10",	"VM_11",	"VM_12",	"VM_2",	"VM_3",	"VM_4",	"VM_5",	"VM_6",	"VM_7", "VM_8", "VM_9"), trim.internal = TRUE, subtree = FALSE,
                 root.edge = 0,collapse.singles = TRUE,
                 interactive = FALSE)
tree_drop$tip.label
# Load the database "state" and specify which states would be pictured
obj<-phylo.to.map(tree_drop,coords,database="state",
                  regions=c("Nevada","California","Oregon"),plot=F)

plot(obj,direction="rightwards",colors=cols,xlim=c(-100,100),ylim=c(-100,100))
plot(obj,type="phylogram",asp=1)
```
