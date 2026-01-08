## Running ANGSD nucleotide diversity

### create a "sites" file to calculate diversity on

```sh
vcftools
--vcf IN.vcf
--kept-sites
--out OUT
```

Need to remove first line that has column headers to index in ANGSD

```sh
sed -i '1d' OUT.kept.sites
```

Then need to index sites file

```sh
startConda
conda activate angsd
angsdPath
```

```sh
~/angsd/angsd sites index OUT.kept.sites
```

Example of running the batch diversity script for ACTH. **RUNNING FROM THE BWA DIRECTORY DUE TO PATHS TO .BAM FILES**

```sh
nohup ../scripts/diversity_parallel_batch.sh wi_test ../angsd/wi_6rare 30 > ../angsd/wi_6rare_batch.log 2>&1 &
```

Reminder that this also works...

```sh
nohup JOBS=30 ../scripts/diversity_parallel_batch.sh wi_test ../angsd/wi_6rare > ../angsd/wi_6rare_batch.log 2>&1 &
```

Where:

+ `wi_test` is a directory within `bwa`that holds each individual population list of rarefied .bam files
    + 1683 total population lists (each with 6 randomly chosen individuals per)
    + 99 random samplings per 17 populations
+ `../angsd/wi_6rare` is the output directory where I'm sending all log files and the final $\pi$ estimates
+ `30` is the number of parallel jobs (number of populations run at once, as the angsd pipeline explicitly limits each population to 1 thread)
