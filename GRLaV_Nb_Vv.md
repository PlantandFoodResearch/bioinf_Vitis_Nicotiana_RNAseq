GRLaV_Nb_Vv
========================================================
author: 
date: 
width: 1920
height: 1080

First Slide
========================================================

For more details on authoring R presentations please visit <https://support.rstudio.com/hc/en-us/articles/200486468>.

- Bullet 1
- Bullet 2
- Bullet 3

Load required libraries
========================================================

```r
library(edgeR)
```

Import count data 
========================================================

Both *N. benthamiana* and *V. vitis* infected with GRLaV, RNAseq performed and read counts
per gene produced. Import the read count `csv` files.


```r
GRLaV3_vs_Nb <- 
  read.delim('/workspace/hradxj/karmun_awesome_experiment/010.edgeR_Nb/GRLaV3_Nb_EdgeR.tab',header=TRUE)

GRLaV3_vs_Vv <- 
  read.delim('/workspace/hradxj/karmun_awesome_experiment/011.edgeR_Vv/GRLaV3_Vv_EdgeR.tab',header=TRUE)
```

Nb expression: data housekeeping
========================================================

Row names are not unique for some reason. Poor annotations.
Assign unique row names using the preexisting names, but making them
unique where required by adding ".1, .2" etc.


```r
rownames(GRLaV3_vs_Nb) <- make.names(GRLaV3_vs_Nb[,1], unique=TRUE)
```

Remove the redundant column where the row names came from


```r
GRLaV3_vs_Nb[,1] <- NULL
```

Nb expression: Replicates and normalisation
========================================================

Set the replicate relationship.. columns 1 and 2 are reps of each other

```r
replicate_relationship <- factor(c(1,1,2))
```
Construct the DGEList object

```r
de_genes <- DGEList(counts=GRLaV3_vs_Nb,group=replicate_relationship,)
```
Calculate normalisation factors

```r
de_genes <- calcNormFactors(de_genes)
```
Examine normalisation factors as a sanity check

```r
de_genes$samples
```

```
                group lib.size norm.factors
Benth.Healthy.1     1 28216160    0.7439606
Benth.Healthy.2     1  2448120    1.1444213
Benth.Infected      2 26897662    1.1745301
```
We can see that Benth.Healthy2 had a lot less sequencing.

Nb expression: Filtering
========================================================
Nonsensical results can happen if we do not remove genes that are not biologically relevant because of low expression. Filter out low expressed genes on the basis of the count per million. This filters out genes that are less than `cpmlimit`, and that are not expressed in `libexplimit` number of libraries.

```r
cpmlimit <- 10
libexplimit <- 3
keep <- rowSums(cpm(de_genes)>cpmlimit) >= libexplimit
```
How many genes were kept?

```r
sum(keep)
```

```
[1] 8982
```
Calculate normalisation factors again and view them

```r
de_genes_kept <- de_genes[keep, , keep.lib.sizes=FALSE]
de_genes_kept <- calcNormFactors(de_genes_kept)
de_genes_kept$samples
```

```
                group lib.size norm.factors
Benth.Healthy.1     1 25563447    0.7319365
Benth.Healthy.2     1  2135625    1.1659241
Benth.Infected      2 23103153    1.1718077
```

Nb expression: Experimental design
========================================================

Set experiment design

```r
design <- model.matrix(~replicate_relationship)
```

Estimate dispersion

```r
de_genes_kept <- estimateDisp(de_genes_kept,design)
```

Nb expression: Test for differential expression I
========================================================

Now we test for differential expression using two methods.
Use the **likelihood ratio test** to test for differentially expressed genes

From the edgeR manual:
>"edgeR offers many variants on analyses. The glm approach is more popular than the classic
>approach as it offers great flexibilities. There are two testing methods under the glm framework:
>likelihood ratio test and quasi-likelihood F-test. The quasi-likelihood method is highly
>recommended for differential expression analyses of bulk RNA-seq data as it gives stricter
>error rate control by accounting for the uncertainty in dispersion estimation. The likelihood
>ratio test can be useful in some special cases such as single cell RNA-seq and datasets with
>no replicates."

Use the glm to obtain DE genes

```r
fitglm <- glmFit(de_genes_kept,design)
lrtglm <- glmLRT(fitglm,coef=2)
topTags(lrtglm, n=10)
```

```
Coefficient:  replicate_relationship2 
                           logFC   logCPM       LR       PValue        FDR
Niben101Scf11283g00002  2.935930 5.681613 21.41761 3.693637e-06 0.03317625
Niben101Scf01731g10001  2.893885 5.854957 19.48912 1.011744e-05 0.04543741
Niben101Scf04083g04043  2.460543 5.274167 17.81218 2.438176e-05 0.07299900
Niben101Scf04209g01002 -3.745765 7.514631 15.82984 6.930075e-05 0.13254312
Niben101Scf04535g02012 -3.627732 8.324320 15.52067 8.160796e-05 0.13254312
Niben101Scf07579g03004  2.247612 5.278974 14.95325 1.102078e-04 0.13254312
Niben101Scf02994g03003  2.503910 5.842941 14.94932 1.104376e-04 0.13254312
Niben101Scf04091g00013 -4.687121 7.885056 14.82355 1.180522e-04 0.13254312
Niben101Scf10735g00018  1.935028 4.851141 13.43418 2.470805e-04 0.24658635
Niben101Scf00449g06008 -3.753175 7.534617 12.52186 4.022185e-04 0.36127265
```

Nb expression: Test for differential expression II
========================================================

Use the **quasi-likelihood model** to obtain differentially expressed genes.
From edgeR manual: 
>"While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
>F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
>provides more robust and reliable error rate control when the number of replicates is small"


```r
fitqlm <- glmQLFit(de_genes_kept, design)
qlftest <- glmQLFTest(fitqlm, coef=2)
topTags(qlftest, n=10)
```

```
Coefficient:  replicate_relationship2 
                           logFC   logCPM        F       PValue       FDR
Niben101Scf11283g00002  2.937170 5.681613 27.42108 0.0008881159 0.9991337
Niben101Scf01731g10001  2.894968 5.854957 25.38431 0.0011253182 0.9991337
Niben101Scf04083g04043  2.460479 5.274167 21.60388 0.0018230061 0.9991337
Niben101Scf04209g01002 -3.745686 7.514631 21.51936 0.0018441156 0.9991337
Niben101Scf04535g02012 -3.627748 8.324320 20.05876 0.0022629039 0.9991337
Niben101Scf02994g03003  2.502932 5.842941 19.42909 0.0024802905 0.9991337
Niben101Scf07579g03004  2.247682 5.278974 18.12272 0.0030223297 0.9991337
Niben101Scf04091g00013 -4.687687 7.885056 16.41105 0.0039811917 0.9991337
Niben101Scf10735g00018  1.934698 4.851141 15.84878 0.0043782090 0.9991337
Niben101Scf06890g00001  2.287987 6.560274 15.27046 0.0048402274 0.9991337
```
