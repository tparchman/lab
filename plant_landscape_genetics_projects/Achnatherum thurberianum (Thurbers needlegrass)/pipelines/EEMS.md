---
title: ""
author: ""
date: ""
output: html_document
---
#### EEMS

```{r eval=FALSE}
####
https://github.com/dipetkov/eems.git
#This repository contains an implementation of the EEMS method for analyzing and visualizing spatial population structure from geo-referenced genetic samples
``` 

Download Eigen 3.2.2 (Eigen does not requires installation)  
Download Boost 1.57 (other version might not be compatibles)  

For Boost installation run the following:
```{r eval=FALSE}
sed -e '1 i#ifndef Q_MOC_RUN' \
    -e '$ a#endif'            \
    -i boost/type_traits/detail/has_binary_operator.hpp &&

./bootstrap.sh --prefix=/usr &&
./b2 stage threading=multi link=shared

sudo su ./b2 install threading=multi link=shared
``` 

After that, go to the Makefile from the runeems_snps/src file and modify EIGEN_INC, BOOST_LIB, and BOOST_INC  

The run, make linux  

EEMS for SNPs: Three input files are required, data.diffs, data.coord, data.outer  

1. data.diffs: the matrix of average pairwise genetic differences. This can be computed with bed2diffs (check before running it the differences among both version of bed2diffs availables)

Go to the bed2diffs/src folder  

bed2diffs uses the libplinkio library to read genotype data stored in plink binary format. To install libplinkio, first clone the GitHub repository and get the latest version (commit 781e9ee37076)  

```{r eval=FALSE}
git clone https://github.com/mfranberg/libplinkio
cd libplinkio
git checkout 781e9ee37076
```

Install libplinkio to a custom location /path/to/plinkio
```{r eval=FALSE}
../configure --prefix=/path/to/plinkio
make && make check && make install
```

Update the PLINKIO location in the Makefile, both in the src and in the src-without-openmp directories  
Compile using make linux  

Usage (in the src folder): 
```{r eval=FALSE}
./bed2diifs_v1 --bfile ./nameofthefiles* --nthreads 2
```

*name of the bed fam and ped files without the extension

Then you got the dissimilarity matrix. Load's in R:

```{r eval=FALSE}
diffs <- read.table("./test/example-SNP-major-mode.diffs")
diffs
#the matrix should be a nonnegative, symmetric with zeros on the diagonal 
```

2. data.coord: the sample coordinates (two coordinates per sample, one sample per line)
3. datapath.outer: the habitat coordinates (as a sequence of vertices that outline a closed polygon). You can get these coordinates in Google Maps API v3 Tool  

Then, EEMS requires a configuration file with "ini" extension, for example:

```{r eval=FALSE}
datapath = ./data/inputdata
mcmcpath = ./data/outputdata
nIndiv = 300
nSites = 3000
nDemes = 200
diploid = false
numMCMCIter = 2000000
numBurnIter = 1000000
numThinIter = 9999

./runeems_snps --params configurationfile.ini --seed 123 (randome seed, it's optional)
```

EEMS results visualization:

```{r eval=FALSE}
#install_github("dipetkov/eems/plotting/rEEMSplots")
#remotes::install_github("dipetkov/reemsplots2")
setwd("../EEMS/")
library("rgdal")   
library(devtools)    
library(rEEMSplots)
library(rworldxtra)
library(reemsplots2)
library(ggplot2)
mcmcpath = "library(rEEMSplots)"
mcmcpath = "../EEMS/"
plotpath = "../EEMS/"
datapath <- file.path("../EEMS/", "filesname")

coord__long_lat <- read.table(paste0('file', ".coord"))
projection_none<-"+proj=longlat +datum=WGS84"
projection_mercator<-"+proj=merc +datum=WGS84"

# Produce the five EEMS figures, with default values for all optional parameters
eems.plots(mcmcpath = "../EEMS/",plotpath = paste0("-default"),longlat = TRUE)

# Add a high-resolution geographic map
eems.plots(
mcmcpath = "../EEMS/",
plotpath = paste0("-geographic-map"),
longlat = TRUE,
projection.in = projection_none,
projection.out = projection_mercator,
add.map = TRUE,
col.map = "black",
lwd.map = 5)

# Add the map explicitly by passing the shape file
map_world <- getMap()
map_na <- map_world[which(map_world@data$continent == "North America"), ]
eems.plots(
mcmcpath = "/home/caro/Documentos/EEMS/eems/runeems_snps/src/condiploidesy300demes/",
plotpath = paste0("-shapefile"),
longlat = TRUE,
m.plot.xy = {plot(map_na, col = NA, add = TRUE)},
q.plot.xy = {plot(map_na, col = NA, add = TRUE)})
```