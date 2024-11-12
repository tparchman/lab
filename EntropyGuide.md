## RUNNING ENTROPY

To activate conda (once installed) on ponderosa

```sh
source ~/anaconda3/bin/activate
```

Create conda env just for entropy (if it doesn't already exist)

```sh
conda create -n entropy
conda activate entropy
```

Either add bioconda channel and install

```sh
conda config --append channels bioconda
conda install popgen-entropy
```

Or install via:

```sh
conda install bioconda::popgen-entropy
```

As of 10/02/2024 this will install gsl 2.7. Need to install gsl 2.6 via (which will give you `libgsl.so.25`):

```sh
conda install gsl=2.6
```

Make an entropy directory within your project directory

```sh
mkdir entropy
cd PATH/entropy
cp KRLA.01.30.2.15.a70.i40.recode.vcf ./
```

This first script would produce a file named `KRLA.01.30.2.15.a70.i40.recode.mpgl` <br>
These scripts can be found on ponderosa in `/mnt/parchman_lab/tfaske/denovo/src/perl_scripts/`

```sh
perl PATH/vcf2mpgl_universal.pl KRLA.01.30.2.15.a70.i40.recode.vcf
```

```sh
perl PATH/gl2genest_universal.pl KRLA.01.30.2.15.a70.i40.recode.mpgl mean
```

Do the following in R:

```r
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(MASS)
library(LEA)

# change as appropriate
setwd('PATH')
```

```r
makePopId <- function(fileIndv){
  PopIDdf = read.table(fileIndv, sep="\t") %>%
    as.data.frame() %>%
    rename(All = V1) %>%
    mutate(Population = str_split_i(All, '_', 2),
           ID = str_split_i(All, '_', 3))
  return(PopIDdf)
}

PCA_entropy <- function(g){
  colmean = apply(g, 2, mean, na.rm = T)
  normalize = matrix(nrow = nrow(g), ncol = ncol(g))
  af = colmean/2
  for (m in 1:length(af)){
    nr = g[,m]-colmean[m]
    dn = sqrt(af[m]*(1-af[m]))
    normalize[,m] = nr/dn
  }
  normalize[is.na(normalize)] = 0
  method1 = prcomp(normalize, scale. = F,center = F)
  pca_df = method1$x[,1:27]
  return(pca_df)
}
```

Need `.indv` output from vcftools. Will come with generating a `.012` file.

```r
PopID <- makePopId('KRLA.01.30.2.15.a70.i40.recode.vcf.012.indv')

g <- read.table('pntest_mean_KRLA.01.30.2.15.a70.i40.recode.txt', header = F)

pca_df <- PCA_entropy(t(g)) %>%
    .[,1:10] %>%
    cbind(PopID)
```

Do kmeans clustering & LDA with desired k-values

```r
writeLDAfile <- function(pcaDF, k){
    kCluster = kmeans(pcaDF[,1:5], k, iter.max = 10, nstart = 10, algorithm = 'Hartigan-Wong')
    ldaOut = lda(x = pcaDF[,1:5], grouping = kCluster$cluster, CV = T)
    write.table(round(ldaOut$posterior, 5),
                file = paste('ldak', as.character(k), '.txt', sep = ''),
                quote = F, row.names = F, col.names = F)
}

writeLDAfile(pca_df, 2)
writeLDAfile(pca_df, 3)
writeLDAfile(pca_df, 4)
writeLDAfile(pca_df, 5)
writeLDAfile(pca_df, 6)
writeLDAfile(pca_df, 7)
writeLDAfile(pca_df, 8)
```

Create header for `.mpgl` file

```r
PopID_list <- paste(PopID$Pop, PopID$ID, sep = '_')

header <- data.frame(dims = NA, PopID_list)

df <- t(header)
dims <- paste(dim(g)[2], dim(g)[1], sep = " ")

df[1,1] <- dims

write.table(df, 'entropy_header.txt',
            sep = " ", na = "",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
```

Back to terminal

```sh
cat entropy_header.txt KRLA.01.30.2.15.a70.i40.recode.mpgl > entropy.mpgl
```
Example for running it on k2... <br>
In Trevor's workflow, he runs 5 chains per k-value and averages. For POMA I ran 3 chains for k 2-7, which is still probably overkill. <br>
<br>
Also, `-r` is only necessary if running different chains of the same k-value at the same time (i.e. you need different seeds to initialize the MCMC). Otherwise the seed is set based on the clock.

```sh
entropy -i entropy.mpgl -o entropy_k2_c1.hdf5 -r RANDOMINT?? -n 2 -l 60000 -b 10000 -t 10 -s 50 -e .01 -k 2 -q ldak2.txt -m 1 -w 0
```

Using python to work with `estpost`

```py
import sys
import pandas as pd
import numpy as np
import scipy as sp
import glob
import re
import random
```

```py
np.set_printoptions(precision=8) # increases float print option
pd.set_option("display.precision", 8)
```

```py
hdf5_files = !find . -name '*hdf5'
hdf5_files = hdf5_files.sort()
hdf5_files
```

```py
estpost = '~/anaconda3/envs/entropy/bin/estpost.entropy'
```

```py
#make DIC
for i in range(0,len(hdf5_files)):
    f = hdf5_files[i]
    k = f.split('_')[1] #set this 
    c = f.split('_')[2].split('.hdf5')[0]
    #print(k,c)
    dic = "DIC_%s_%s.txt" % (k,c)
    !$estpost $f -s 3 -p deviance > $dic
```

```py
dic_files = !find . -name 'DIC*'
dic_files
```

```py
for d in dic_files:
    !cat $d
    print('\n')
```

```py
dic_list = []
for d in dic_files:
    k = d.split('_k')[1].split('_')[0] #set this 
    c = d.split('_c')[1].split('.txt')[0]
    #print(k,c)
    
    dic = !grep 'DIC' $d
    dic = float(re.search('(\d+.\d+)',str(dic)).group(0))
    #print(dic)
    
    dic_list.append([k,dic,c])
dic_df = pd.DataFrame(dic_list,columns=['k','DIC','chain'])
dic_df.head()
```

```py
dic_df.to_csv('dic_list.csv')
```

```py
dic_sum = dic_df.groupby('k').describe().DIC
```

```py
dic_sum.sort_values('mean')
```

```py
dic_sum.sort_values('mean')
```

```py
dic_sum.to_csv('dic_sum.csv')
```

Generating ancestry coefficients (combining all chains if there are multiple?)

```sh
# ancestry coeffecients 
!$estpost *k2*.hdf5 -p q -s 0 -o q2.txt

!$estpost *k3*.hdf5 -p q -s 0 -o q3.txt

!$estpost *k4*.hdf5 -p q -s 0 -o q4.txt

!$estpost *k5*.hdf5 -p q -s 0 -o q5.txt

!$estpost *k6*.hdf5 -p q -s 0 -o q6.txt
```

```sh
#MCMC diagnostics
!$estpost *k2*.hdf5 -p q -s 4 -o MCMC_k2.txt

!$estpost *k3*.hdf5 -p q -s 4 -o MCMC_k3.txt

!$estpost *k4*.hdf5 -p q -s 4 -o MCMC_k4.txt

!$estpost *k5*.hdf5 -p q -s 4 -o MCMC_k5.txt

!$estpost *k6*.hdf5 -p q -s 4 -o MCMC_k6.txt
```

Generate genotype probability file for k2 (also combining chains if there are multiple?)

```sh
estpost *k2*.hdf5 -p gprob -s 0 -o gprob2.txt
```

Can also generate a gprob file that's an average of gprobs across all run k-values

```py
hdf5_files = []
num_k = [2,3,4,5,6]
num_c = 4
for k in num_k:
    for c in range(1,num_c+1):
        f = '../entropy_k' + str(k) + '_c' + str(c) + '.hdf5'
        hdf5_files.append(f)
hdf5_files
```

```py
gprob_cmd = estpost + ' ' + ' '.join(hdf5_files) + ' -p gprob -s 0 -o ../gprobAll.txt'
```

```sh
gprob_cmd
```