#### map_ggplot.R
#### Mapping and PCA, modified from Faske maps
#### 

#devtools::install_github("dkahle/ggmap")
library(ggmap)
library(devtools)
require(ggplot2)
require(ggsci)
require(ggrepel)
library(data.table)
library(patchwork)

##### Read in map data #####
setwd("/home/caro/Escritorio/")
coord <- read.csv('PHHA_map.csv', header=T)
names(coord) <- c('Name','Pop','Lat','Long')

#################### read in genotype data##########

map <- get_stamenmap(bbox = c(left = -121.1, bottom = 38.0, right = -108.5, top = 49.5),
                     zoom=8,scale = 3,maptype = 'terrain-background')
##zoom controls resolution, 10 takes forever

#view the map to make sure the dimensions are good
ggmap(map) + geom_point(data = coord, aes(x = Long, y = Lat),size=3,pch=21,fill="white",col="black") +
  xlab("Longitude") + ylab("Latitude") +
  coord_map(xlim=c(-121.1,-108.5),ylim=c(38.0,49.5)) #selects the range for x and y

##### Map with labels by Population ####
map <- ggmap(map) +
  geom_point(data = coord, aes(x = Long, y = Lat),size=2,col="black",fill='white',pch=21) +
  geom_label_repel(data = coord, aes(x = Long, y = Lat,label = Pop,fill=Pop),
                   colour = "black", fontface = "bold") +
  scale_fill_d3(palette = 'category20') +
  #scale_fill_manual(name='Population:',values=listofcolors) +
  #above comment changes the color of the fill in the labels with custom colors
  coord_map(xlim=c(-121.1,-108.5),ylim=c(38.0,49.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size=13),
        axis.title = element_text(size = 15, colour="black",face = "bold", vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
map
