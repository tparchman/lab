
### Stacks

```{r eval=FALSE}
module load stacks/1.46
module load gcc/4.8.5  

denovo_map.pl --samples /archive/parchman_lab/rawdata_to_backup/GSAF_11_20_bc_parse/ACTH/ -o /working/mascaro/ACTH/ -T 5 -O /home/working/ACTH/pop_map -m 3 -M 2 -n 2 -S -b 1
```

-samples: Path to the directory of samples (samples will be read from population map)  
-T: Number of threads  
-O: Population map path  
-m: Specify a minimum number of identical reads required to create a stack  
-M: Specify the number of mismatches allowed between loci when processing a single individual (default 2)  
-n: Specify the number of mismatches allowed between loci when building the catalog (default 1)  
-S: Disable recording SQL data in the database  
-b: Batch ID representing this dataset (an integer, e.g. 1, 2, 3)  

pop_map: Population map file, name of the samples (without the file extension), and a number or string to identified the population, separated by a tab (make sure that there are not hidden/invisible characters in the file)  
To make sure of that, run: 
```{r eval=FALSE}
head population_map od -c
```

The file should looks like:

AT_AH_01 | AH  
AT_AH_02 | AH  
AT_AH_03 | AH  
AT_AH_04 | AH  

To get a vcf file output, run:

```{r eval=FALSE}
populations -b 1 -P /working/mascaro/ACTH -t 12 -M /working/mascaro/ACTH/pop_map --max_obs_het 0.65 -r 0.80 --vcf
```

-max-obs-het: specify a maximum observed heterozygosity required to process a nucleotide site at a locus  
-r: minimum percentage of individuals in a population required to process a locus for that population  

### vcf statistics with vcftools

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
##### Calculate site quality
```{r eval=FALSE}
vcftools --gzvcf VCF --site-quality --out $OUT
```
##### Calculate proportion of missing data per individual
```{r eval=FALSE}
vcftools --gzvcf VCF --missing-indv --out $OUT
```
##### Calculate proportion of missing data per site
```{r eval=FALSE}
vcftools --gzvcf VCF --missing-site --out $OUT
```

### Examining statistics in R
```{r eval=FALSE}
library(tidyverse)

##Phred encoded site quality
var_qual <- read_delim("site-quality.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

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

##Heterozygosity and inbreeding coefficient per individual
ind_het <- read_delim("heterozigosity.txt.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
```

### vcf filter with vcftools
```{r eval=FALSE}
vcftools --vcf batch_sin_.vcf --maf 0.03 --max-missing 0.8 --minQ 10 --min-meanDP 5 --max-meanDP 16 --minDP 5 --maxDP 16 --recode --out variantes_filtro.vcf
```

-maf: minor-allele frequency thresholds  
-max-missing: set minimum missing data (0 is totally missing, 1 is none missing)  
-minQ: minimum quality score required for a site to pass our filtering threshold  
-min-meanDP: the minimum mean depth for a site  
-max-meanDP: the maximum mean depth for a site (good rule: the mean depth x2)  
-minDP: the minimum depth allowed for a genotype  
-maxDP: the maximum depth allowed for a genotype  
