Basic filtering steps (how I typically do it)

Basic rules:  
1. You always need an input argument:
    + `--gzvcf inputFile.vcf.gz` if your input is compressed
    + `--vcf inputFile.vcf` if uncompressed

2. You need to identify an output format (choose 1):
    + If outputting to a new VCF file, use both of these so you retain information like depth
        + `--recode` and `--recode-INFO-all`  
        output file suffix = ***.recode.vcf**
    + Various "summary" files
        + `--missing-indv` (proportion of total sites missing for each individual)  
        output file suffix = ***.imiss**
        + `--depth` (average depth for each individual)  
        output file suffix = ***.idepth**
        + `--freq2` (allelic frequencies for each site)  
        output file suffix = ***.frq**
        + `--site-mean-depth`  
        output file suffix = ***.ldepth.mean**

3. You need to name the output file
    + via `--out outputFile`
    You do **NOT** need the suffix (will be added based on output format)

4. **IF** you're filtering (actually removing sites or individuals), you also need:
    + `--remove-filtered-all`  
    Not needed if you're simply getting a summary file

I typically name things based on any filtering paramters that have been used. For instance, a fully filtered file might be named something like **KRLA.bithin.q30.i70.a80.maf03.maxdp15** which means the following arguments would have been used:  
+ **KRLA** - 4-letter code for the taxa
+ **bithin** - biallelic only, thinned to 100, indels removed
+ **q30** - site quality score > 30
+ **i70** - individuals removed if they are missing >70% of sites
+ **a80** - only sites that are present in >=80% of individuals
+ **maf03** - only sites wth minor allele frequency >=0.03
+ **maxdp15**  - only sites with mean depth values (over included individuals) <=15

I typically take the following steps:  
1. Filter for biallelic SNPs only, thinned, and call quality >30. Produce a downsized VCF.
2. Summarize individual missingness
3. Produce downsized VCF based on removed individuals
4. Summarize site depth
5. Produce downsized VCF based on max depth
6. Summarize SNP totals across combinations of potential site-missingness and MAF thresholds
7. Try some site-missingness/MAF filtering values and assess PCA
8. **OPTIONAL:** filter on additional parameters (maybe necessary on case by case basis)

Initially filter on 3 things:

+ SNPs only (no indels) `--remove-indels`
+ biallelic sites only `--min-alleles 2` and `--max-alleles 2`
+ 1 SNP per contig (helps account for LD, our reads are typically >100 BP length) `--thin 100`
+ **OPTIONAL (but I typically do)** - filter on call quality. Our call qualities are typically very good so it doesn't get rid of much but might as well cut out the worst stuff `--minQ 30`


```sh
vcftools \
--gzvcf KRLA.vcf.gz \
--remove-indels \
--min-alleles 2 \
--max-alleles 2 \
--thin 100 \
--minQ 30 \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out KRLA.bithin.q30
```
Produces file called **KRLA.bithin.q30.recode.vcf**

Make a summary file to look at individual missingness

```sh
vcftools \
--vcf KRLA.bithin.q30.recode.vcf \
--missing-indv \
--out KRLA.bithin.q30
```
Produces file called **KRLA.bithin.q30.imiss**


Transfer to local system (do this from the local folder you want to move **INTO**)
```sh
scp romero@ponderosa.biology.unr.edu:PATH_TO_SERVER_FOLDER/KRLA.bithin.q30.imiss ./
```

**NEXT STEPS DONE IN R**

These should be all the packages you need to do R stuff in this document:

```r
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(data.table)
library(ggplot2)
library(colorspace)
library(ggrepel)
library(patchwork)
```

Make some color palette with enough colors for the number of populations you have. Here's 37 to play around with.
```r
col37 <- c("#a1def0", "#b32728", "#4ed31b", "#bf209e", "#3f862d", 
           "#b32df9", "#99d683", "#74398b", "#ead624", "#4853fc", 
           "#fe8f06", "#395f97", "#1cf1a3", "#ed8bc7", "#7c440e", 
           "#c6c0fe", "#444a1e", "#efbba2", "#369094", "#fa756b",
           "#88e99a", "#256b33", "#6dc5dd", "#154975", "#628df2",
           "#874fbf", "#e84fe1", "#eec8f1", "#752e4f", "#c7799a",
           "#2a2bf0", "#aebf8a", "#743502", "#bce333", "#e72525",
           "#fea53b", "#ce1365")
```
This will make a figure where you can look at how different populations vary in terms of missing data. Identify some potential thresholds of missingness where you might want to filter out individuals.

```r
ind_miss  <- read_delim("KRLA.bithin.q30.imiss", 
                        delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), 
                        skip = 1) %>%
    # next line should work for most sample naming conventions, but adjust as needed
    mutate(Pop = str_split_i(ind, '_', 2))

imissFig <- ggplot(ind_miss, 
                   aes(x = fmiss, y = Pop, fill = Pop, color = Pop)) +
  geom_jitter(pch = 21, size = 2,
              position = position_jitter(height = 0.2)) +
  scale_fill_manual(values = col37) +
  scale_color_manual(values = darken) +
  theme_bw() +
  theme(legend.position = 'none')

imissFig
```
**NOW BACK TO COMMAND LINE**

Use this command to create a list of individuals to filter out. In this case you're identifying all individuals missing more than 70% of the total sites. **MAKE SURE YOU ALSO CHANGE THE NAME OF OUTPUT .txt FILE FOR ORGANIZATION**

```sh
awk '$5 > 0.7 {print $1}' KRLA.imiss | tail -n +2 > indmiss70.txt
```
Now filter out those identified individuals and make a downsized VCF via the following:

```sh
vcftools \
--vcf KRLA.bithin.q30.recode.vcf \
--remove indmiss70.txt
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out KRLA.bithin.q30.i70
```
Next, I like to look at depth across the individuals that have reasonable coverage. Create an output summary file via:

```sh
vcftools \
--vcf KRLA.bithin.q30.i70.recode.vcf \
--site-mean-depth \
--out KRLA.bithin.q30
```
Produces file called **KRLA.bithin.q30.i70.ldepth.mean**

Move this output file to local system via `scp` as shown previously, then **RUN THE FOLLOWING IN R**

```r
var_depth <- read_delim("KRLA.bithin.q30.i70.ldepth.mean", delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

vdepthFig <- ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  # for ddradseq, average depth is going to be <20X, so this will remove exteme high depth outliers
  theme_light() + xlim(0, 30) +
  ggtitle('Variant Depth (cutoff above 30)') +
  theme(plot.title = element_text(hjust = 0.5))

vdepthFig
```

Also good to look at distribution of depth via:

```r
summary var_depth$mean_depth
```

There's more that goes into this, but I typically will do a max depth cutoff ~15 if the mean depth is under 7. If you have higher depth, doing someting like 2X your mean is fine. There's other more specific formulas out in the literature.

**BACK TO COMMAND LINE**

Filter out site above a maximum depth threshold and create a downsized VCF via:

```sh
vcftools \
--vcf KRLA.bithin.q30.i70.recode.vcf \
--max-meanDP 15
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out KRLA.bithin.q30.i70.maxdp15
```

Next create a summary file of allelic frequencies by site. This will allow us to look how much SNP retention we can expect across a range of filtering parameters for site-missingness and MAF.

```sh
vcftools \
--vcf KRLA.bithin.q30.i70.maxdp15.recode.vcf \
--freq2 \
--out KRLA.bithin.q30.i70.maxdp15
```

Produces a file called **KRLA.bithin.q30.i70.maxdp15.frq**

Again, move to local system via `scp` and do the following in R

Code below essentially calculates the number of loci that would be **retained** for a given combination of site-missingness and MAF. The figure produced will allow you to understand how much data will be lost/gained by choosing more inclusive/restrictive filters. I generally try to maximize loci while simultaneously minimize missing data. The MAF is more dependent on the number of total samples and number of individuals per population. Typically I'll choose something in the 0.02 - 0.05 range.

```r
var_freq <- read_delim("KRLA.bithin.q30.i70.maxdp15.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_freq$maf <- var_freq %>% dplyr::select(a1, a2) %>% apply(1, function(z) min(z))
var_freq <- mutate(var_freq, fmiss = (max(nchr)-nchr)/max(nchr))

maf_set <- matrix(nrow = 51, ncol = 20)
maf_set[,1] <- seq(0, 0.50, 0.01)
amiss_set <- seq(95, 5, -5)
for(i in 1:nrow(maf_set)){
  for(j in 1:length(amiss_set)){
    maf_set[i,(j+1)] <- as.integer(nrow(filter(var_freq, maf > ((i-1)/100) & fmiss < (amiss_set[j]/100))))
  }
}
ind_with <- as.data.frame(maf_set)

fNames <- vector()
for(i in 1:length(amiss_set)){
  fNames[i] = paste('fmiss', as.character(amiss_set[i]), sep = '_')
}
colnames(ind_with)[2:20] <- fNames

comb <- ind_with %>%
  pivot_longer(cols = 2:20, names_to = 'fmiss', values_to = 'snps')

param_compare <- ggplot(comb, aes(x = V1, y = snps, group = fmiss, label = fmiss, color = fmiss)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 0.5, 0.02)) +
  scale_y_continuous(breaks = seq(0, 100000, 10000)) +
  coord_cartesian(xlim = c(0, 0.13),
                  ylim = c(0, 60000)) +
  geom_label_repel(size = ifelse(comb$V1 == 0.01, 5, 0),
                   alpha = ifelse(comb$V1 == 0.01, 1, 0),
                   color = 'black',
                   max.overlaps = 40,
                   label.padding = 0.1, box.padding = 0.65) +
  scale_color_moma_d(palette_name = 'Warhol') +
  labs(x = 'MAF',
       y = 'Total Loci') +
  ggtitle('Loci Retention Rate') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')

param_compare
```

Decide on a combination of MAF/missigness values and return to command line. Note that to filter on missingness, you need to use the inverse of the numbers shown in that above figure. For example, `--max-missing 0.8` means that you are only keeping sites that have data in **at least 80%** of all individuals (this is the same as saying that missingness **cannot exceed 20%**).

From this point, I typically just output 012 files so that I can explore what the filtered SNPs look like in PCA. I'll do one final `--recode` output only once I've decided on final filtering parameters.

To create these 012 files and test parameters, do something like the following:

```sh
vcftools \
--vcf KRLA.bithin.q30.i70.maxdp15.recode.vcf \
--max-missing 0.8
--maf 0.05
--remove-filtered-all \
--012 \
--out KRLA.bithin.q30.i70.maxdp15.a80.maf05
```

This will produce three files with the suffixes `.012`, `.012.indv`, and `.012.pos`. You can move all of these at once using escaped wildcards like so:

```sh
scp romero@ponderosa.biology.unr.edu:PATH_TO_SERVER_FOLDER/KRLA.bithin.q30.i70.maxdp15.a80.maf05\*012\* ./
```

Finally, you can explore the filtered SNPs in PCA by using a combination of the following functions.

First, a function to format your population ID file:

```r
makePopId <- function(fileIndv){
  PopIDdf = read.table(fileIndv, sep="\t") %>%
    as.data.frame() %>%
    rename(All = V1) %>%
    # note that this will work for the most basic sample ID formats but will need to be modified for projects that use other sample ID patterns...
    mutate(Population = str_split_i(All, '_', 2))
  return(PopIDdf)
}
```

Next, a function that will generate two PCAs side by side (PC1xPC2 and PC3xPC4).

```r
pca1234FUN <- function(file012, popID, colorSet, grouping){
  df012 = fread(file012, sep = '\t', data.table = F)[,-1] %>%
    as.matrix() %>%
    apply(2, function(d) gsub(-1, NA, d, fixed=TRUE)) %>%
    apply(2, function(d) as.numeric(d))
  colmean = apply(df012, 2, mean, na.rm=TRUE)
  normalize = matrix(nrow = nrow(df012), ncol = ncol(df012))
  af = colmean/2
  for (m in 1:length(af)){
    nr = df012[ ,m]-colmean[m]
    dn = sqrt(af[m]*(1-af[m]))
    normalize[ ,m] = nr/dn
  }
  normalize[is.na(normalize)] = 0
  pca012 = prcomp(normalize, scale. = FALSE, center = FALSE)
  pca012_loading = pca012$x
  pc1var = formatC(((summary(pca012)$importance[2,1])*100), digits = 1, format = 'f')
  pc2var = formatC(((summary(pca012)$importance[2,2])*100), digits = 1, format = 'f')
  pc3var = formatC(((summary(pca012)$importance[2,3])*100), digits = 1, format = 'f')
  pc4var = formatC(((summary(pca012)$importance[2,4])*100), digits = 1, format = 'f')
  pc5var = formatC(((summary(pca012)$importance[2,5])*100), digits = 1, format = 'f')
  pc6var = formatC(((summary(pca012)$importance[2,6])*100), digits = 1, format = 'f')
  pca_df = pca012_loading[,1:6]
  pca_df = cbind(popID, pca_df) 
  pca12_pop = ggplot(data = pca_df, aes(x = PC1, y = PC2, 
                                        fill = .data[[grouping]],
                                        color = .data[[grouping]])) +
    geom_point(pch = 21, size = 3) + 
    xlab(paste("PC",1," (",pc1var,"%)",sep="")) +
    ylab(paste("PC",2," (",pc2var,"%)",sep="")) +
    scale_fill_manual(name = paste(grouping, ':', sep = ''), values = colorSet) +
    scale_color_manual(name = paste(grouping, ':', sep = ''), values = darken(colorSet, 0.3)) +
    # scale_fill_moma_c('Ernst') +
    theme_bw() + 
    theme(legend.position = 'none',
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16, color = "black",
                                    face = "bold", vjust = 1),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  pca34_pop = ggplot(data = pca_df, aes(x = PC3, y = PC4,
                                        fill = .data[[grouping]],
                                        color = .data[[grouping]])) +
    geom_point(pch = 21, size = 3) + 
    xlab(paste("PC",3," (",pc3var,"%)",sep="")) +
    ylab(paste("PC",4," (",pc4var,"%)",sep="")) +
    scale_fill_manual(name = paste(grouping, ':', sep = ''), values = colorSet) +
    scale_color_manual(name = paste(grouping, ':', sep = ''), values = darken(colorSet, 0.3)) +
    theme_bw() + 
    theme(legend.position = 'none',
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16, color = "black",
                                    face = "bold", vjust = 1),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  finalPatch = pca12_pop + pca34_pop + plot_layout(nrow = 1)
  return(finalPatch)
}
```

Finally, a function that will make a fireworks plat (with labels for each population). You can designate which axes you want to look at for this.

```r
fwStd <- function(file012, popID, xAxis, yAxis, colorSet, grouping){
  df012 = fread(file012, sep = '\t', data.table = F)[,-1] %>%
    as.matrix() %>%
    apply(2, function(d) gsub(-1, NA, d, fixed=TRUE)) %>%
    apply(2, function(d) as.numeric(d))
  colmean = apply(df012, 2, mean, na.rm=TRUE)
  normalize = matrix(nrow = nrow(df012), ncol = ncol(df012))
  af = colmean/2
  for (m in 1:length(af)){
    nr = df012[ ,m]-colmean[m]
    dn = sqrt(af[m]*(1-af[m]))
    normalize[ ,m] = nr/dn
  }
  normalize[is.na(normalize)] = 0
  pca012 = prcomp(normalize, scale. = FALSE, center = FALSE)
  pcXvar = formatC(((summary(pca012)$importance[2,as.integer(str_sub(xAxis, 3, 3))])*100), 
                   digits = 1, format = 'f')
  pcYvar = formatC(((summary(pca012)$importance[2,as.integer(str_sub(yAxis, 3, 3))])*100), 
                   digits = 1, format = 'f')
  pca012_loading = pca012$x
  pca_df = pca012_loading[,1:6]
  pca_df = cbind(popID, pca_df)
  pca_mean <- pca_df %>% 
    group_by(.data[[grouping]]) %>% 
    summarize(PC1_mean=mean(PC1),PC2_mean=mean(PC2),PC3_mean=mean(PC3),
              PC4_mean=mean(PC4),PC5_mean=mean(PC5),PC6_mean=mean(PC6))
  pca_df <- left_join(pca_df,pca_mean)
  pca_fw <- ggplot(data = pca_mean, 
                   aes(x = .data[[paste(xAxis,'mean',sep='_')]],
                       y = .data[[paste(yAxis,'mean',sep='_')]],
                       fill = .data[[grouping]],
                       color = .data[[grouping]], 
                       label = .data[[grouping]])) +
    geom_segment(data=pca_df,
                 aes(x = .data[[xAxis]],y = .data[[yAxis]],
                     xend = .data[[paste(xAxis,'mean',sep='_')]],
                     yend = .data[[paste(yAxis,'mean',sep='_')]]),
                 linewidth = 0.3) +
    geom_label_repel(color = 'black',
                     max.overlaps = 40, label.padding = 0.1, box.padding = 0.65,
                     label.size = 0.1,
                     alpha = 0.5,
                     size = 1,
                     segment.size = 0.1) +
    geom_point(size = 4, 
               stroke = 0.5, pch = 21) + 
    xlab(paste("PC",(str_sub(xAxis, 3, 3))," (",pcXvar,"%)",sep="")) +
    ylab(paste("PC",(str_sub(yAxis, 3, 3))," (",pcYvar,"%)",sep="")) +
    scale_fill_manual(name = paste(grouping, ':', sep = ''), values = colorSet) +
    scale_color_manual(name = paste(grouping, ':', sep = ''), values = darken(colorSet, 0.3)) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text = element_text(size = 10),
          axis.title.y.left = element_text(size = 16,
                                           face = 'bold',
                                           vjust = 0),
          axis.title.x.bottom = element_text(size = 16,
                                             face = 'bold',
                                             vjust = 1),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          legend.text = element_text(size = 13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = 'pt'))
  return(pca_fw)
}
```

Once these functions exist in your environment, you can easily re-create PCA/fireworks plots across different filtering parameters via the following:

```r
# first one takes the .012.indv file
Apop <- makePopId('KRLA.bithin.q30.i70.maxdp15.a80.maf05.012.indv')
Apca <- pca1234FUN('KRLA.bithin.q30.i70.maxdp15.a80.maf05.012', Mpop, col37, 'Population')
# can make arguments explicit
Afw12 <- fwStd(file012 = 'KRLA.bithin.q30.i70.maxdp15.a80.maf05.012',
               popID = Apop, 
               xAxis = 'PC1', 
               yAxis = 'PC2',
               grouping = 'Population', 
               colorSet = col37)
# or just list in correct order
Afw34 <- fwStd('KRLA.bithin.q30.i70.maxdp15.a80.maf05.012',
               Apop, 
               'PC3'
               'PC4',
               'Population', 
               col37)


A_4panel <- Apca + Afw12 + Afw34 + plot_layout(nrow = 2)

A_4panel
```


# Below are old filtering notes from Trevor if useful...



### Filtering on Individuals

* **Coverage:** Also can be thought of as depth. See 3mapping. Calculated on bam files. Average read count per locus per individual.
* **Missing:** Proportion of missing data allowed across all loci for individual. Common and high in GBS/RADseq data. Kinda an issue all around. Many methods, including PCA (all ordination methods), require a complete matrix with no missing data. Additionally, PCA will cluster by missing data with individuals with higher missing data clustering closer to the center and get this "fan" effect. Can be the same for coverage too. This (among other reasons) is why people use a variance-covariance matrix of genetic data to do ordinations. Other methods involve imputation. This can be fancy and use phased haplotype data OR simply, when you z-score, (g - mean(g))/sd(g), your genotype data across each locus, you make all missing data equal to 0 or Mean (i.e., the global allele frequency). There's more to this standardization, see [Patterson et al. 2006](https://dx.plos.org/10.1371/journal.pgen.0020190) for more info. See PCAsim_ex in examples directory for showing all these issues. This is another reason to use entropy. Entropy is a hierarchical bayesian model so it gets an updated genotype estimate for each missing value based on genotype likelihoods across loci, individuals, and the allele frequency of the cluster/deme that individual assigns to.

### Filtering on Loci

* **Biallelic:** Only keep biallelic SNPs. Multiallelic SNPs are rare at the time scale we work (Citation??) and also, mathematical nightmare and we have enough data so just ignore. Everyone does unless deep time phylogenetics.
* **thin:** Keeps one locus within a specified range. Not 100% how it decides with one to keep. I think it's on quality or depth. This is a necessary step as loci in close physical are prone to sequencing error and linkage disequalibrium (LD) confounds many different population genetic parameters. For de novo reference assemblies, we thin to 100 as contigs/reads are ~92 bp in length. This keeps one locus per contig to control for LD and sequencing error, which is really common in pop gen and necessary for many analyses.
* **max-missing** = max proportion of missing data per locus
* **MAF** = minor allele frequency. Proportion of individuals a alternate allele needs to be present in order for that locus to be kept as a SNP. (e.g. maf = 0.02 for 250 individuals means that an alternate allele needs to be present in at least 5 individuals to be kept) Many papers have shown this is a big issue in clustering and demography (Citation). We do this a second time near the end if we removed individuals during missing data filtering.
* **Mean Depth:** Average allelic depth or read depth per locus. Too low could be sequencing error, too high could be PCR replication artifact (Citation).
* **Qual:** Locus Quality. You can look up the math. Usually above 20-30 is good but given our coverage and number of individuals, we can usually go way higher.
* **Fis:** Inbreeding coefficient. This is a contentous topic. This has to do with paralogs or paralogous loci. This is where loci map to multiple regions of the genome. Issues in highly repeative genomes. Usually leads to an excess of heterozygotes. Filtering on negative Fis can help. See these two McKinney papers [1](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12763), and [2](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12613). Katie and others in the lab use his package called HDPlot to deal with this.