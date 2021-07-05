### Haplotype networks:

1.Convert vcf format to fasta using vcf2phylip.py (https://github.com/edgardomortiz/vcf2phylip)

```{r eval=FALSE}
python vcf2phylip.py --i acth.vcf --fasta
```
 2.Trimming the alignment with trimAl (http://trimal.cgenomics.org/_media/tutorial.pdf):
  
```{r eval=FALSE}
./trimal -in acth.fasta -gappyout -out acth_trimmed.fasta
```

3.Haplotype network in R (ape, pegas):

```{r eval=FALSE}
input <- "acth_trimmed.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
h <- sort(h, what = "label")
(net <- pegas::haploNet(h))
ind.hap<-with(
stack(setNames(attr(h, "index"), rownames(h))),
table(hap=ind, pop=rownames(d)[values]))

plot(net, size=attr(net, "freq"), scale.ratio=0.2, pie=ind.hap)
legend(-8, 0, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)

```
