### Admixture

#### Generating the input file

```{r eval=FALSE}
./plink2 --vcf $FILE.vcf --allow-extra-chr --out fileset  

./plink2 --vcf $FILE.vcf --allow-extra-chr --maf 0.05  --make-bed --out fileset  
```

ADMIXTURE does not accept chromosome names that are not human chromosomes. Change the first column by 0 for the bim, bed, and fam files:

```{r eval=FALSE}
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim
```

#### Running Admixture
```{r eval=FALSE}
for i in {2..18}
do
 ./admixture --cv FILE.bed $i > log${i}.out
done
```

##### Get the cv errors to identify the best value of k clusters which is the value with lowest cross-validation error  

```{r eval=FALSE}
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > FILE.cv.error
```

##### Plot the results using admixture.R script from https://github.com/speciationgenomics/scripts  

```{r eval=FALSE}
wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r
chmod +x plotADMIXTURE.r

Rscript admixture.R -p filesnames -i individual.populations.txt -k 8 -m 8 -l AH,AS,BM,BV,DC,DH,EW,FR,GB,HO,JC,LV,MD,PL,PT,SC,SS,VM
```

-p filesnames without the extension  
-i individual populations file (individuals in a column, pop in other one)  
-k maximum number of K's  
-m minimum number of K's  
-l population name  


#### Evaluate the results of the admixture analysis with EvalMix: http://www.popgen.dk/software/index.php/EvalAdmix

#### EvalAdmix installation  

```{r eval=FALSE}
git clone https://github.com/GenisGE/evalAdmix.git
cd evalAdmix
make

./evalAdmix -plink files -fname file.P -qname file.Q -P 20 
```
-plink: Genotypes file prefix files    
-fname: Ancestral Populations frequency file file.P  
-qname: Admixture proportions file file.Q
-P: Computer cores 

#### Plot the results in R
```{r eval=FALSE}
source("visFuns.R")

# Read population labels and estimated admixture proportions
pop<-read.table("file.fam")
q<-read.table("file.Q")

# Order according to population and plot the ADMIXTURE reults
ord<-orderInds(pop = as.vector(pop[,2]), q = q)
barplot(t(q)[,ord],col=2:10,space=0,border=NA,xlab="Individuals",ylab="Demo2 Admixture proportions for K=3")
text(tapply(1:nrow(pop),pop[ord,2],mean),-0.05,unique(pop[ord,2]),xpd=T)
abline(v=cumsum(sapply(unique(pop[ord,2]),function(x){sum(pop[ord,2]==x)})),col=1,lwd=1.)

r<-as.matrix(read.table("output.corres.txt"))

# Plot correlation of residuals
plotCorRes(cor_mat = r, pop = as.vector(pop[,2]), ord=ord, title="Evaluation of 1000G admixture proportions with K=3", max_z=0.1, min_z=-0.1)
```
