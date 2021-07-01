### ACTH 

Cleaned and demultiplexed .fastq files by individual:

```{r eval=FALSE}
/archive/parchman_lab/rawdata_to_backup/GSAF_11_20_bc_parse/ACTH
```

### Calling variants using dDocent (http://www.ddocent.com/): main stepts 

```{r eval=FALSE}
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n ddocent_env ddocent
source activate ddocent_env
dDocent

/working/mascaro/acth/

mkdir reads

scp /archive/parchman_lab/rawdata_to_backup/GSAF_11_20_bc_parse/ACTH/*fastq.gz /working/mascaro/acth/reads

#dDocent requires files named *F.fq.gz
conda install -c bioconda rename
rename 's/fastq/F.fq/g' *fastq.gz
```

### Optimize denovo reference: 

#### 1. Place a subset of individuals of the total data set (5 individuals per populations) in refOpt.

#### 2. Run ReferenceOpt.sh: This script assembles references across cutoff values and then maps X random samples and evaluates mappings to the reference, along with number of contigs and coverage.


```{r eval=FALSE}
mkdir refOpt

./ReferenceOpt.sh 4 10 4 10 0.9 16 SE

#This would loop across cutoffs of 4-10 using a similarity of 90% for clustering, parellized across 16 processors, using SE assembly technique.
# The output is stored in a file called mapping.results

```
#### 3. Visualize data in kopt.data: Plot values for each k1,k2 combination across similarity thresholds and pick a similarity threshold at the point of inflection on the curve:

```{r eval=FALSE}
library(ggplot2)

data.table <- read.table("kopt.data", header = FALSE, col.names= c("k1","k2","Similarity", "Contigs"))
data.table$K1K2 <- paste(data.table$k1, data.table$k2, sep=",")
df=data.frame(data.table)
df$K1K2 <- as.factor(df$K1K2)

p <- ggplot(df, aes(x=Similarity, y=Contigs, group=K1K2)) + scale_x_continuous(breaks=seq(0.8,0.98,0.01)) + geom_line(aes(colour = K1K2))
p

# Similarity threshold of 0.94 

```

### Run RefMapOpt.sh using 0.94 as similarity threshold:

Reads need to be trimmed before run RefMapOpt!

Pick optimal k1,k2 cutoffs (Ideally, you want to maximize properly paired mappings and coverage while minimizing mismatched reads.)

### Run dDocent on this subset with the correct assembly parameters (skipping mapping and snp calling.)

Copy the reference.fasta file from this RefOpt directory to your main working directory.

### Run dDocent on your full data set, skipping trimming and assembly.

