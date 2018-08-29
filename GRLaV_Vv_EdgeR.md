RNA-Seq aligned against V. vinifera
================
Marcus Davy, Dan Jones

<!-- Note: Rendering of urls in github rmarkdown will not work -->
> The master version of this document is on *[github](https://github.com/PlantandFoodResearch/bioinf_Vitis_Nicotiana_RNAseq)*. This printing is shasum version *[c04984e](https://github.com/PlantandFoodResearch/bioinf_Vitis_Nicotiana_RNAseq/commit/c04984e)* from that repository.

> Rmd =&gt; R files rendered with;

    ## knitr::purl("GRLaV_Vv_EdgeR.Rmd", documentation=2L)

Count data
----------

Checking names in dataset do not contain duplicated names, or missing names. A missing annotation could imply an off by *n* allocation of gene names to count data, which would make all gene look ups after the first position of the mislabeled record meaningless.

``` r
filename <- '/workspace/hradxj/karmun_awesome_experiment/011.edgeR_Vv/GRLaV3_Vv_EdgeR.tab'

GRLaV3_vs_Vv <- read.delim(filename, header=TRUE)
dim(GRLaV3_vs_Vv)
```

    ## [1] 19759     4

``` r
#GRLaV3_vs_Vv2 <- read.table(filename, header=TRUE, sep="\t")

## Sanity check names
## Check for missing Gene names
GRLaV3_vs_Vv[which(GRLaV3_vs_Vv$Gene %in% ""), ]
```

    ## [1] Gene             Grape.Healthy    Grape.Infected.1 Grape.Infected.2
    ## <0 rows> (or 0-length row.names)

``` r
## Duplicated names
badnames <- names(table(GRLaV3_vs_Vv$Gene))[table(GRLaV3_vs_Vv$Gene)>1]
GRLaV3_vs_Vv[which(GRLaV3_vs_Vv$Gene %in% badnames), ]
```

    ## [1] Gene             Grape.Healthy    Grape.Infected.1 Grape.Infected.2
    ## <0 rows> (or 0-length row.names)

``` r
# Row names are not unique for some reason. Crap annotations again.
# Assign unique row names using the preexisting names, but making them
# unique where required by adding ".1, .2" etc.
rownames(GRLaV3_vs_Vv) <- make.names(GRLaV3_vs_Vv[,1], unique=TRUE)

## Sanity check
GRLaV3_vs_Vv[which(GRLaV3_vs_Vv$Gene %in% badnames),]
```

    ## [1] Gene             Grape.Healthy    Grape.Infected.1 Grape.Infected.2
    ## <0 rows> (or 0-length row.names)

``` r
# Remove the redundant column where the row names came from
GRLaV3_vs_Vv[,1] <- NULL
```

Data exploration
----------------

``` r
par(pty="s")

suppressWarnings(pairs(log2(GRLaV3_vs_Vv), labels=colnames(GRLaV3_vs_Vv), upper.panel=panel.points, lower.panel=panel.corr,
      pch=16, cex=0.4, asp=1))
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/exploration-1.png)

``` r
suppressWarnings(pairs(log2(GRLaV3_vs_Vv), labels=colnames(GRLaV3_vs_Vv), upper.panel=panel.points, lower.panel=panel.corr,
      maTransform=TRUE, xlim=c(0,20), ylim=c(-8,8), pch=16, cex=0.4, asp=1))
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/exploration-2.png)

``` r
#plotMDS(log2(GRLaV3_vs_Vv+0.5), top=1000)


## How many samples are 0 in All samples?
# Filtering these has no effect on cpm/rpkm calculations
apply(GRLaV3_vs_Vv==0, 1, all)  %>%  table
```

    ## .
    ## FALSE  TRUE 
    ## 16986  2773

``` r
# set the replicate relationship.. columns 1 and 2 are reps of each other
## Healthy vs Infected
replicate_relationship <- factor(c(1,2,2))
# construct the DGEList object
all_genes <- DGEList(counts=GRLaV3_vs_Vv, group=replicate_relationship)

for( i in seq(ncol(all_genes))) {
  plotMD.DGEList(all_genes, column=i)
  abline(h=0, col="grey", lty=2)
}
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/exploration-3.png)![](GRLaV_Vv_EdgeR_files/figure-markdown_github/exploration-4.png)![](GRLaV_Vv_EdgeR_files/figure-markdown_github/exploration-5.png)

``` r
breaks <- seq(0, ceiling(max(log2(all_genes$counts))))
tweak  <- 1

## Distribution of counts
for(i in colnames(all_genes)) {
  hist(log2(all_genes$counts+tweak)[, i], breaks=breaks, main=i, xlab="log2(counts)")
}
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/exploration-6.png)![](GRLaV_Vv_EdgeR_files/figure-markdown_github/exploration-7.png)![](GRLaV_Vv_EdgeR_files/figure-markdown_github/exploration-8.png)

`Grape.Infected.1` replicate is different from the other samples (including the other biological replicate)

The spike in the histograms are the zero counts in each sample.

Filtering
---------

Library normalization factors estimated using trimmed mean of M values (TMM).

Alternative less aggressive filtering of only genes where all data was zero counts was used. In this experiment over filtering has the additional side effect of effecting the mixture distribution of modeled raw p-values which should be distributed Ho: unif(0,1) , Ha: exp(rate).

``` r
# Calculate normalisation factors
all_genes <- calcNormFactors(all_genes, method="TMM")

# Examine normalisation factors as a sanity check
all_genes$samples
```

    ##                  group lib.size norm.factors
    ## Grape.Healthy        1  5845892    1.0267447
    ## Grape.Infected.1     2 18046557    0.9686044
    ## Grape.Infected.2     2  1521236    1.0055208

``` r
# We can see that Benth.Healthy, Grape.Infected.2 had a lot less sequencing.

# Nonsensical results can happen if we do not remove genes
# that are not biologically relevant because of low expression
# Filter out low expressed genes on the basis of the count per million
# This filters out genes that are less than cpmlimit, and 
# that are not expressed in libexplimit number of libraries
cpmlimit    <- 10
libexplimit <- 2
keep <- rowSums(cpm(all_genes)>cpmlimit) >= libexplimit

## Visualize
hist(log2(rowSums(cpm(all_genes))))
abline(v=log2(cpmlimit), col="red", lty=2)
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/filtering-1.png)

``` r
## Sanity check how many are kept
table(keep)
```

    ## keep
    ## FALSE  TRUE 
    ## 10380  9379

``` r
par(pty="s")

## Visualize filtering
suppressWarnings(pairs(log2(GRLaV3_vs_Vv), labels=colnames(GRLaV3_vs_Vv), upper.panel=panel.points, lower.panel=panel.corr, pch=16, cex=0.1, asp=1, col=c("black","grey")[(keep)+1], main="Filtering effect"))
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/filtering-2.png)

``` r
suppressWarnings(pairs(log2(GRLaV3_vs_Vv), labels=colnames(GRLaV3_vs_Vv), upper.panel=panel.points, lower.panel=panel.corr, maTransform=TRUE, xlim=c(0,20), ylim=c(-8,8), pch=16, cex=0.1, asp=1, col=c("black", "grey")[(keep)+1], main="MAplot Filtering effect"))
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/filtering-3.png)

``` r
## Distribution of counts after filtering
for(i in colnames(all_genes)) {
  hist(log2(all_genes$counts+tweak)[keep , i], breaks=breaks, main=i, xlab="log2(counts)")
}
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/filtering-4.png)![](GRLaV_Vv_EdgeR_files/figure-markdown_github/filtering-5.png)![](GRLaV_Vv_EdgeR_files/figure-markdown_github/filtering-6.png)

``` r
## Alternative filter which overrides line:89
# keep <- !apply(GRLaV3_vs_Vv==0, 1, all)

## Sanity check how many are kept
# table(keep)

all_genes_kept <- all_genes[which(keep), keep.lib.sizes=FALSE]

## Sanity check
dim(all_genes)
```

    ## [1] 19759     3

``` r
dim(all_genes_kept)
```

    ## [1] 9379    3

``` r
# Calculate normalisation factors again
all_genes_kept <- calcNormFactors(all_genes_kept, method="TMM")

## Sanity check
apply(GRLaV3_vs_Vv,2, sum)
```

    ##    Grape.Healthy Grape.Infected.1 Grape.Infected.2 
    ##          5845892         18046557          1521236

``` r
calcNormFactors(all_genes, method="TMM")
```

    ## An object of class "DGEList"
    ## $counts
    ##                   Grape.Healthy Grape.Infected.1 Grape.Infected.2
    ## GSVIVG01012261001            14               37                1
    ## GSVIVG01012259001             6               18                1
    ## GSVIVG01012257001           311             1151              151
    ## GSVIVG01012255001           755             1819              234
    ## GSVIVG01012253001             1                0                0
    ## 19754 more rows ...
    ## 
    ## $samples
    ##                  group lib.size norm.factors
    ## Grape.Healthy        1  5845892    1.0267447
    ## Grape.Infected.1     2 18046557    0.9686044
    ## Grape.Infected.2     2  1521236    1.0055208

``` r
all_genes_kept$samples
```

    ##                  group lib.size norm.factors
    ## Grape.Healthy        1  5628229    1.0290413
    ## Grape.Infected.1     2 17591761    0.9782519
    ## Grape.Infected.2     2  1480593    0.9933825

``` r
## Grape.Infected.2 is an order of magnitude smaller!
apply(GRLaV3_vs_Vv,2, sum)[2] / apply(GRLaV3_vs_Vv,2, sum)[3]
```

    ## Grape.Infected.1 
    ##         11.86309

Of note here is the `Grape.Infected.1` sample has 11.8630883 times the library size of `Grape.Infected.2`. Why are the two libraries so different?

Filtering out all gene models where no count is observed has no effect on library size estimation.

edgeR models
------------

``` r
# Todo: Need to look at underlying bam mapping scores...


# set experiment design
design <- model.matrix(~replicate_relationship)
# estimate dispersion
all_genes_kept <- estimateDisp(all_genes_kept, design)


# Now we test for differential expression using two methods.
# From the edgeR manual:
# "edgeR offers many variants on analyses. The glm approach is more popular than the classic
# approach as it offers great flexibilities. There are two testing methods under the glm framework:
# likelihood ratio test and quasi-likelihood F-test. The quasi-likelihood method is highly
# recommended for differential expression analyses of bulk RNA-seq data as it gives stricter
# error rate control by accounting for the uncertainty in dispersion estimation. The likelihood
# ratio test can be useful in some special cases such as single cell RNA-seq and datasets with
# no replicates."

# Use the glm to obtain DE genes
fitglm <- glmFit(all_genes_kept, design)
lrtglm <- glmLRT(fitglm, coef=2)
topN   <- topTags(lrtglm, n=nrow(lrtglm))


summary(decideTests(lrtglm))
```

    ##        replicate_relationship2
    ## Down                       222
    ## NotSig                    8885
    ## Up                         272

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/modelling-1.png)

``` r
openDevice("png", file.path(imgDir, "maplot_lrt"))
```

    ## [1] "images/Vv/maplot_lrt.png"

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
dev.off()
```

    ## png 
    ##   2

``` r
openDevice("tiff", file.path(imgDir, "maplot_lrt"))
```

    ## [1] "images/Vv/maplot_lrt.tif"

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
dev.off()
```

    ## png 
    ##   2

``` r
threshold <- 0.05
de_genes  <- topN$table[topN$table$FDR<threshold,]

## DE genes
nrow(de_genes)
```

    ## [1] 494

``` r
## First 10
head(de_genes, n=10)
```

    ##                       logFC    logCPM        LR       PValue          FDR
    ## GSVIVG01028856001 -7.155443  9.571178 111.47956 4.645693e-26 4.357196e-22
    ## GSVIVG01017960001 -6.620759  8.084244  81.51203 1.741978e-19 8.169006e-16
    ## GSVIVG01035900001 -4.764818  8.128116  61.55656 4.302057e-15 1.344966e-11
    ## GSVIVG01026014001 -4.723862  6.672261  59.16300 1.451376e-14 3.403113e-11
    ## GSVIVG01027028001 -4.637934 10.682202  51.96331 5.654749e-13 1.060718e-09
    ## GSVIVG01033388001 -4.384228  6.434971  51.27756 8.018683e-13 1.253454e-09
    ## GSVIVG01017087001 -4.532545  6.127439  50.19270 1.393654e-12 1.867297e-09
    ## GSVIVG01019407001 -4.770217 11.925204  48.92392 2.660855e-12 3.119520e-09
    ## GSVIVG01011496001 -4.188728  8.576402  47.39930 5.790325e-12 6.034162e-09
    ## GSVIVG01005149001 -3.985182  6.472926  43.41146 4.435836e-11 4.160371e-08

``` r
## estimate pi0
pvals <- topN$table$PValue

hist(pvals)
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/modelling-2.png)

``` r
## Investigating abnormal right hand peak in pvalue histogram
# ind <- which(pvals>=1)
# rownames(topN$table)[ind]
# GRLaV3_vs_Vv[rownames(topN$table)[ind],]

if(est_pi0) {
  pi0_hat <- limma::convest(pvals)
  print(pi0_hat)
}
```

    ## [1] 0.8772344

``` r
# Use the quasi-likelihood model to obtain DE genes.
# From edgeR manual: 
#"While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small"
fitglm <- glmQLFit(all_genes_kept, design)
qlfglm <- glmQLFTest(fitglm, coef=2)
topN2  <- topTags(qlfglm, n=nrow(qlfglm))


summary(decideTests(qlfglm))
```

    ##        replicate_relationship2
    ## Down                       135
    ## NotSig                    9105
    ## Up                         139

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/modelling-3.png)

``` r
openDevice("png", file.path(imgDir, "maplot_ql"))
```

    ## [1] "images/Vv/maplot_ql.png"

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
dev.off()
```

    ## png 
    ##   2

``` r
openDevice("tiff", file.path(imgDir, "maplot_ql"))
```

    ## [1] "images/Vv/maplot_ql.tif"

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
dev.off()
```

    ## png 
    ##   2

``` r
threshold <- 0.05
de_genes2 <- topN2$table[topN2$table$FDR<threshold,]

## DE genes
nrow(de_genes2)
```

    ## [1] 274

``` r
## First 10
head(de_genes2, n=10)
```

    ##                       logFC   logCPM         F       PValue          FDR
    ## GSVIVG01028856001 -7.154303 9.571178 136.37533 1.334937e-08 0.0001252038
    ## GSVIVG01017960001 -6.644303 8.084244  96.01189 1.208468e-07 0.0005667113
    ## GSVIVG01035900001 -4.764070 8.128116  81.90847 3.185049e-07 0.0008096809
    ## GSVIVG01026014001 -4.724850 6.672261  80.81827 3.453165e-07 0.0008096809
    ## GSVIVG01033388001 -4.384892 6.434971  67.56378 1.002384e-06 0.0018802711
    ## GSVIVG01017087001 -4.533882 6.127439  63.41298 1.451461e-06 0.0022688752
    ## GSVIVG01011496001 -4.188637 8.576402  60.42399 1.919098e-06 0.0025526102
    ## GSVIVG01026942001 -3.844252 7.354024  58.21496 2.376654e-06 0.0025526102
    ## GSVIVG01005149001 -3.984124 6.472926  57.15654 2.639316e-06 0.0025526102
    ## GSVIVG01004352001 -3.799065 7.575296  56.84949 2.721623e-06 0.0025526102

``` r
## estimate pi0
pvals2 <- topN$table$PValue

hist(pvals2)
```

![](GRLaV_Vv_EdgeR_files/figure-markdown_github/modelling-4.png)

``` r
## Investigating abnormal right hand peak in pvalue histogram
# ind <- which(pvals>=1)
# rownames(topN$table)[ind]
# GRLaV3_vs_Vv[rownames(topN$table)[ind],]

if(est_pi0) {
  pi0_hat2 <- limma::convest(pvals2)
  print(pi0_hat2)
}
```

    ## [1] 0.8772344
