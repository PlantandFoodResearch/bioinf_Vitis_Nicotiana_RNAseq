RNA-Seq aligned against N. benthamiana
================
Marcus Davy, Dan Jones

<!-- Note: Rendering of urls in github rmarkdown will not work -->
> The master version of this document is on *[github](https://github.com/PlantandFoodResearch/bioinf_Vitis_Nicotiana_RNAseq)*. This printing is shasum version *[c04984e](https://github.com/PlantandFoodResearch/bioinf_Vitis_Nicotiana_RNAseq/commit/c04984e)* from that repository.

> Rmd =&gt; R files rendered with;

    ## knitr::purl("GRLaV_Nb_EdgeR.Rmd", documentation=2L)

Count data
----------

Checking names in dataset do not contain duplicated names, or missing names. A missing annotation could imply an off by *n* allocation of gene names to count data, which would make all gene look ups after the first position of the mislabeled record meaningless.

``` r
filename <- '/workspace/hradxj/karmun_awesome_experiment/010.edgeR_Nb/GRLaV3_Nb_EdgeR.tab'

GRLaV3_vs_Nb <- read.delim(filename, header=TRUE)
dim(GRLaV3_vs_Nb)
```

    ## [1] 59814     4

``` r
#GRLaV3_vs_Nb2 <- read.table(filename, header=TRUE, sep="\t")

## Sanity check names
## Check for missing Gene names
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% ""), ]
```

    ## [1] Gene            Benth.Healthy.1 Benth.Healthy.2 Benth.Infected 
    ## <0 rows> (or 0-length row.names)

``` r
## Duplicated names
badnames <- names(table(GRLaV3_vs_Nb$Gene))[table(GRLaV3_vs_Nb$Gene)>1]
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% badnames), ]
```

    ##                         Gene Benth.Healthy.1 Benth.Healthy.2
    ## 10644 Niben101Scf01064g06002               0               0
    ## 10645 Niben101Scf01064g06002              26               0
    ##       Benth.Infected
    ## 10644              0
    ## 10645             48

``` r
# Row names are not unique for some reason. Crap annotations again.
# Assign unique row names using the preexisting names, but making them
# unique where required by adding ".1, .2" etc.
rownames(GRLaV3_vs_Nb) <- make.names(GRLaV3_vs_Nb[,1], unique=TRUE)

## Sanity check
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% badnames),]
```

    ##                                            Gene Benth.Healthy.1
    ## Niben101Scf01064g06002   Niben101Scf01064g06002               0
    ## Niben101Scf01064g06002.1 Niben101Scf01064g06002              26
    ##                          Benth.Healthy.2 Benth.Infected
    ## Niben101Scf01064g06002                 0              0
    ## Niben101Scf01064g06002.1               0             48

``` r
# Remove the redundant column where the row names came from
GRLaV3_vs_Nb[,1] <- NULL
```

Data exploration
----------------

``` r
par(pty="s")

suppressWarnings(pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr,
      pch=16, cex=0.4, asp=1))
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/exploration-1.png)

``` r
suppressWarnings(pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr,
      maTransform=TRUE, xlim=c(0,20), ylim=c(-8,8), pch=16, cex=0.4, asp=1))
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/exploration-2.png)

``` r
#plotMDS(log2(GRLaV3_vs_Nb+0.5), top=1000)


## How many samples are 0 in All samples?
# Filtering these has no effect on cpm/rpkm calculations
apply(GRLaV3_vs_Nb==0, 1, all)  %>%  table
```

    ## .
    ## FALSE  TRUE 
    ## 44511 15303

``` r
# set the replicate relationship.. columns 1 and 2 are reps of each other
## Healthy vs Infected
replicate_relationship <- factor(c(1,1,2))
# construct the DGEList object
all_genes <- DGEList(counts=GRLaV3_vs_Nb, group=replicate_relationship)

for( i in seq(ncol(all_genes))) {
  plotMD.DGEList(all_genes, column=i)
  abline(h=0, col="grey", lty=2)
}
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/exploration-3.png)![](GRLaV_Nb_EdgeR_files/figure-markdown_github/exploration-4.png)![](GRLaV_Nb_EdgeR_files/figure-markdown_github/exploration-5.png)

``` r
breaks <- seq(0, ceiling(max(log2(all_genes$counts))))
tweak  <- 1

## Distribution of counts
for(i in colnames(all_genes)) {
  hist(log2(all_genes$counts+tweak)[, i], breaks=breaks, main=i, xlab="log2(counts)")
}
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/exploration-6.png)![](GRLaV_Nb_EdgeR_files/figure-markdown_github/exploration-7.png)![](GRLaV_Nb_EdgeR_files/figure-markdown_github/exploration-8.png)

`Healthy.2` replicate is different from the other samples (including the other biological replicate)

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

    ##                 group lib.size norm.factors
    ## Benth.Healthy.1     1 28216160    0.7439606
    ## Benth.Healthy.2     1  2448120    1.1444213
    ## Benth.Infected      2 26897662    1.1745301

``` r
# We can see that Benth.Healthy2 had a lot less sequencing.

# Nonsensical results can happen if we do not remove genes
# that are not biologically relevant because of low expression
# Filter out low expressed genes on the basis of the count per million
# This filters out genes that are less than cpmlimit, and 
# that are not expressed in libexplimit number of libraries
cpmlimit    <- 1
libexplimit <- 1
keep <- rowSums(cpm(all_genes)>cpmlimit) >= libexplimit

## Visualize
hist(log2(rowSums(cpm(all_genes))))
abline(v=log2(cpmlimit), col="red", lty=2)
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/filtering-1.png)

``` r
## Sanity check how many are kept
table(keep)
```

    ## keep
    ## FALSE  TRUE 
    ## 26920 32894

``` r
par(pty="s")

## Visualize filtering
suppressWarnings(pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr, pch=16, cex=0.1, asp=1, col=c("black","grey")[(keep)+1], main="Filtering effect"))
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/filtering-2.png)

``` r
suppressWarnings(pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr, maTransform=TRUE, xlim=c(0,20), ylim=c(-8,8), pch=16, cex=0.1, asp=1, col=c("black", "grey")[(keep)+1], main="MAplot Filtering effect"))
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/filtering-3.png)

``` r
## Alternative filter which overrides line:89
keep <- !apply(GRLaV3_vs_Nb==0, 1, all)

## Sanity check how many are kept
table(keep)
```

    ## keep
    ## FALSE  TRUE 
    ## 15303 44511

``` r
## Visualize filtering
suppressWarnings(pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr, pch=16, cex=0.1, asp=1, col=c("black","grey")[(keep)+1], main="Filtering effect"))
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/filtering-4.png)

``` r
suppressWarnings(pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr, maTransform=TRUE, xlim=c(0,20), ylim=c(-8,8), pch=16, cex=0.1, asp=1, col=c("black", "grey")[(keep)+1], main="MAplot Filtering effect"))
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/filtering-5.png)

``` r
## Distribution of counts after filtering
for(i in colnames(all_genes)) {
  hist(log2(all_genes$counts+tweak)[keep , i], breaks=breaks, main=i, xlab="log2(counts)")
}
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/filtering-6.png)![](GRLaV_Nb_EdgeR_files/figure-markdown_github/filtering-7.png)![](GRLaV_Nb_EdgeR_files/figure-markdown_github/filtering-8.png)

``` r
all_genes_kept <- all_genes[which(keep), keep.lib.sizes=FALSE]

## Sanity check
dim(all_genes)
```

    ## [1] 59814     3

``` r
dim(all_genes_kept)
```

    ## [1] 44511     3

``` r
# Calculate normalisation factors again
all_genes_kept <- calcNormFactors(all_genes_kept, method="TMM")

## Sanity check
apply(GRLaV3_vs_Nb,2, sum)
```

    ## Benth.Healthy.1 Benth.Healthy.2  Benth.Infected 
    ##        28216160         2448120        26897662

``` r
calcNormFactors(all_genes, method="TMM")
```

    ## An object of class "DGEList"
    ## $counts
    ##                        Benth.Healthy.1 Benth.Healthy.2 Benth.Infected
    ## Niben101Ctg00054g00001              55              23             79
    ## Niben101Ctg00074g00004               0               0              0
    ## Niben101Ctg00116g00002              65              17            133
    ## Niben101Ctg00129g00001               0               0              0
    ## Niben101Ctg00141g00002               1               0              1
    ## 59809 more rows ...
    ## 
    ## $samples
    ##                 group lib.size norm.factors
    ## Benth.Healthy.1     1 28216160    0.7439606
    ## Benth.Healthy.2     1  2448120    1.1444213
    ## Benth.Infected      2 26897662    1.1745301

``` r
all_genes_kept$samples
```

    ##                 group lib.size norm.factors
    ## Benth.Healthy.1     1 28216160    0.7439606
    ## Benth.Healthy.2     1  2448120    1.1444213
    ## Benth.Infected      2 26897662    1.1745301

``` r
## Benth.Healthy.2 is an order of magnitude smaller!
apply(GRLaV3_vs_Nb,2, sum)[1] / apply(GRLaV3_vs_Nb,2, sum)[2]
```

    ## Benth.Healthy.1 
    ##        11.52564

Of note here is the `Benth.Healthy.1` sample has 11.5256442 times the library size of `Benth.Healthy.2`. Why are the two libraries so different?

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
    ## Down                       107
    ## NotSig                   44355
    ## Up                          49

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/modelling-1.png)

``` r
openDevice("png", file.path(imgDir, "maplot_lrt"))
```

    ## [1] "images/Nb/maplot_lrt.png"

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
dev.off()
```

    ## png 
    ##   2

``` r
openDevice("tiff", file.path(imgDir, "maplot_lrt"))
```

    ## [1] "images/Nb/maplot_lrt.tif"

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

    ## [1] 156

``` r
## First 10
head(de_genes, n=10)
```

    ##                            logFC   logCPM       LR       PValue
    ## Niben101Scf00107g03008  6.871699 4.682465 63.44407 1.649850e-15
    ## Niben101Scf05044g02012  6.622191 2.769140 43.62414 3.979007e-11
    ## Niben101Scf00837g08001 -6.892553 4.177794 39.84457 2.749970e-10
    ## Niben101Scf03937g00009  5.188343 3.759200 39.48756 3.301559e-10
    ## Niben101Scf06977g00014 -6.624423 6.394987 35.68045 2.324855e-09
    ## Niben101Scf03169g00010 -5.974855 4.342707 34.04949 5.372795e-09
    ## Niben101Scf01475g00019 -6.064870 3.542031 32.66200 1.096600e-08
    ## Niben101Scf03506g03001  4.587488 3.224802 32.41610 1.244513e-08
    ## Niben101Scf02971g01006 -9.719592 1.913830 31.37670 2.125159e-08
    ## Niben101Scf08683g00001 -5.387247 4.306713 29.82370 4.731716e-08
    ##                                 FDR
    ## Niben101Scf00107g03008 7.343648e-11
    ## Niben101Scf05044g02012 8.855480e-07
    ## Niben101Scf00837g08001 3.673892e-06
    ## Niben101Scf03937g00009 3.673892e-06
    ## Niben101Scf06977g00014 2.069632e-05
    ## Niben101Scf03169g00010 3.985808e-05
    ## Niben101Scf01475g00019 6.924314e-05
    ## Niben101Scf03506g03001 6.924314e-05
    ## Niben101Scf02971g01006 1.051033e-04
    ## Niben101Scf08683g00001 2.106134e-04

``` r
## estimate pi0
pvals <- topN$table$PValue

hist(pvals)
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/modelling-2.png)

``` r
## Investigating abnormal right hand peak in pvalue histogram
# ind <- which(pvals>=1)
# rownames(topN$table)[ind]
# GRLaV3_vs_Nb[rownames(topN$table)[ind],]

if(est_pi0) {
  pi0_hat <- limma::convest(pvals)
  print(pi0_hat)
}
```

    ## [1] 0.9872359

``` r
## Obtaining genes that show zero differential expression
high_expressed_no_de_genes <- subset(topN$table, logFC < 0.01 & logFC > -0.01 & logCPM > 6)
dim(high_expressed_no_de_genes)
```

    ## [1] 28  5

``` r
knitr::kable(high_expressed_no_de_genes)
```

|                        |       logFC|     logCPM|         LR|     PValue|        FDR|
|------------------------|-----------:|----------:|----------:|----------:|----------:|
| Niben101Scf01196g00003 |   0.0096889|   6.043966|  0.0001645|  0.9897653|  0.9995753|
| Niben101Scf02972g01003 |   0.0087496|   8.294383|  0.0001248|  0.9910879|  0.9995753|
| Niben101Scf02164g04001 |   0.0081702|   7.644829|  0.0001102|  0.9916243|  0.9995753|
| Niben101Scf09908g00005 |  -0.0080815|   7.401203|  0.0001082|  0.9916993|  0.9995753|
| Niben101Scf03169g05002 |   0.0079898|   9.577414|  0.0001015|  0.9919604|  0.9995753|
| Niben101Scf35531g00006 |  -0.0076203|   6.287457|  0.0001006|  0.9919990|  0.9995753|
| Niben101Scf02078g05005 |  -0.0070210|   7.569653|  0.0000814|  0.9928001|  0.9995753|
| Niben101Scf10189g02009 |  -0.0068661|   6.773991|  0.0000797|  0.9928784|  0.9995753|
| Niben101Scf02613g03001 |   0.0067699|   8.472678|  0.0000744|  0.9931190|  0.9995753|
| Niben101Scf06516g04018 |   0.0063622|   6.222272|  0.0000703|  0.9933083|  0.9995753|
| Niben101Scf07824g00009 |   0.0058622|   8.631338|  0.0000556|  0.9940511|  0.9996981|
| Niben101Scf00550g00014 |  -0.0050760|   7.052778|  0.0000430|  0.9947692|  0.9996981|
| Niben101Scf01905g00001 |   0.0048720|   7.212177|  0.0000395|  0.9949885|  0.9996981|
| Niben101Scf00063g13005 |   0.0045982|   6.108572|  0.0000369|  0.9951502|  0.9996981|
| Niben101Scf01834g01011 |  -0.0046794|   7.211773|  0.0000364|  0.9951837|  0.9996981|
| Niben101Scf03069g04013 |   0.0043431|   6.330008|  0.0000326|  0.9954447|  0.9997823|
| Niben101Scf01148g02008 |  -0.0040399|   6.015278|  0.0000286|  0.9957313|  0.9997823|
| Niben101Scf02031g00029 |   0.0040084|   6.386766|  0.0000277|  0.9958008|  0.9997823|
| Niben101Scf06579g00007 |  -0.0037393|   6.053308|  0.0000245|  0.9960515|  0.9997823|
| Niben101Scf06964g00007 |   0.0037818|  13.219614|  0.0000221|  0.9962513|  0.9997823|
| Niben101Scf01382g05008 |   0.0034531|   6.927004|  0.0000200|  0.9964299|  0.9997823|
| Niben101Scf01466g01006 |  -0.0023004|   7.136597|  0.0000088|  0.9976312|  0.9998550|
| Niben101Scf09928g01011 |   0.0020133|   7.082736|  0.0000067|  0.9979280|  0.9998741|
| Niben101Scf04703g03003 |   0.0009214|   6.252873|  0.0000015|  0.9990317|  0.9998741|
| Niben101Scf08618g00005 |  -0.0008697|   6.548008|  0.0000013|  0.9990935|  0.9998741|
| Niben101Scf00568g04071 |  -0.0002956|   9.971503|  0.0000001|  0.9997033|  0.9999488|
| Niben101Scf04099g08003 |  -0.0000879|   8.674557|  0.0000000|  0.9999108|  0.9999782|
| Niben101Scf18308g01001 |   0.0000631|   6.922904|  0.0000000|  0.9999347|  0.9999793|

``` r
# Use the quasi-likelihood model to obtain DE genes.
# From edgeR manual: 
#"While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small"
fitqlm  <- glmQLFit(all_genes_kept, design)
qlftest <- glmQLFTest(fitqlm, coef=2)
topN2   <- topTags(qlftest, n=nrow(lrtglm))


summary(decideTests(qlftest))
```

    ##        replicate_relationship2
    ## Down                       138
    ## NotSig                   44319
    ## Up                          54

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/modelling-3.png)

``` r
openDevice("png", file.path(imgDir, "maplot_ql"))
```

    ## [1] "images/Nb/maplot_ql.png"

``` r
plotMD(lrtglm, cex=0.5, main="Infected vs Healthy")
dev.off()
```

    ## png 
    ##   2

``` r
openDevice("tiff", file.path(imgDir, "maplot_ql"))
```

    ## [1] "images/Nb/maplot_ql.tif"

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

    ## [1] 192

``` r
## First 10
head(de_genes2, n=10)
```

    ##                            logFC   logCPM        F       PValue
    ## Niben101Scf00107g03008  6.871688 4.682465 79.73219 4.447203e-19
    ## Niben101Scf00837g08001 -6.892875 4.177794 50.13170 1.458892e-12
    ## Niben101Scf05044g02012  6.622181 2.769140 49.77721 1.747415e-12
    ## Niben101Scf03937g00009  5.188513 3.759200 48.58255 3.210717e-12
    ## Niben101Scf03169g00010 -5.975193 4.342707 43.04653 5.403794e-11
    ## Niben101Scf06977g00014 -6.624507 6.394987 42.22785 8.209053e-11
    ## Niben101Scf01475g00019 -6.065571 3.542031 40.11691 2.414876e-10
    ## Niben101Scf03506g03001  4.587501 3.224802 38.53894 5.414464e-10
    ## Niben101Scf08683g00001 -5.387564 4.306713 37.65900 8.496614e-10
    ## Niben101Scf02175g06018 -5.251927 4.355896 36.51320 1.528372e-09
    ##                                 FDR
    ## Niben101Scf00107g03008 1.979495e-14
    ## Niben101Scf00837g08001 2.592639e-08
    ## Niben101Scf05044g02012 2.592639e-08
    ## Niben101Scf03937g00009 3.572805e-08
    ## Niben101Scf03169g00010 4.810566e-07
    ## Niben101Scf06977g00014 6.089886e-07
    ## Niben101Scf01475g00019 1.535551e-06
    ## Niben101Scf03506g03001 3.012540e-06
    ## Niben101Scf08683g00001 4.202142e-06
    ## Niben101Scf02175g06018 6.802937e-06

``` r
## estimate pi0
pvals2 <- topN2$table$PValue

hist(pvals2)
```

![](GRLaV_Nb_EdgeR_files/figure-markdown_github/modelling-4.png)

``` r
## Investigating abnormal right hand peak in pvalue histogram
# ind <- which(pvals>=1)
# rownames(topN$table)[ind]
# GRLaV3_vs_Nb[rownames(topN$table)[ind],]

if(est_pi0) {
  pi0_hat2 <- limma::convest(pvals2)
  print(pi0_hat2)
}
```

    ## [1] 0.9922079

``` r
## Obtaining genes that show zero differential expression
high_expressed_no_de_genes <- subset(topN2$table, logFC < 0.01 & logFC > -0.01 & logCPM > 6)
dim(high_expressed_no_de_genes)
```

    ## [1] 28  5

``` r
knitr::kable(high_expressed_no_de_genes)
```

|                        |       logFC|     logCPM|          F|     PValue|        FDR|
|------------------------|-----------:|----------:|----------:|----------:|----------:|
| Niben101Scf01196g00003 |   0.0096877|   6.043966|  0.0001977|  0.9887823|  0.9995796|
| Niben101Scf35531g00006 |  -0.0076193|   6.287457|  0.0001188|  0.9913041|  0.9995796|
| Niben101Scf02972g01003 |   0.0087497|   8.294383|  0.0001170|  0.9913710|  0.9995796|
| Niben101Scf09908g00005 |  -0.0080814|   7.401203|  0.0001147|  0.9914568|  0.9995796|
| Niben101Scf02164g04001 |   0.0081699|   7.644829|  0.0001132|  0.9915097|  0.9995796|
| Niben101Scf10189g02009 |  -0.0068664|   6.773991|  0.0000903|  0.9924177|  0.9995796|
| Niben101Scf02078g05005 |  -0.0070206|   7.569653|  0.0000845|  0.9926663|  0.9995796|
| Niben101Scf06516g04018 |   0.0063616|   6.222272|  0.0000835|  0.9927100|  0.9995796|
| Niben101Scf03169g05002 |   0.0079899|   9.577414|  0.0000762|  0.9930367|  0.9995796|
| Niben101Scf02613g03001 |   0.0067697|   8.472678|  0.0000679|  0.9934275|  0.9995807|
| Niben101Scf07824g00009 |   0.0058620|   8.631338|  0.0000494|  0.9943916|  0.9996934|
| Niben101Scf00550g00014 |  -0.0050766|   7.052778|  0.0000474|  0.9945078|  0.9996934|
| Niben101Scf00063g13005 |   0.0045977|   6.108572|  0.0000442|  0.9946957|  0.9996934|
| Niben101Scf01905g00001 |   0.0048738|   7.212177|  0.0000428|  0.9947802|  0.9996934|
| Niben101Scf01834g01011 |  -0.0046792|   7.211773|  0.0000395|  0.9949871|  0.9996934|
| Niben101Scf03069g04013 |   0.0043423|   6.330008|  0.0000384|  0.9950575|  0.9996934|
| Niben101Scf01148g02008 |  -0.0040408|   6.015278|  0.0000345|  0.9953156|  0.9996934|
| Niben101Scf02031g00029 |   0.0040084|   6.386766|  0.0000325|  0.9954533|  0.9996934|
| Niben101Scf06579g00007 |  -0.0037392|   6.053308|  0.0000294|  0.9956733|  0.9996934|
| Niben101Scf01382g05008 |   0.0034532|   6.927004|  0.0000224|  0.9962275|  0.9996934|
| Niben101Scf01466g01006 |  -0.0023010|   7.136597|  0.0000096|  0.9975237|  0.9998394|
| Niben101Scf09928g01011 |   0.0019893|   7.082736|  0.0000073|  0.9978496|  0.9998394|
| Niben101Scf06964g00007 |   0.0037818|  13.219614|  0.0000070|  0.9978935|  0.9998394|
| Niben101Scf04703g03003 |   0.0009220|   6.252873|  0.0000017|  0.9989455|  0.9998702|
| Niben101Scf08618g00005 |  -0.0008699|   6.548008|  0.0000015|  0.9990246|  0.9998702|
| Niben101Scf00568g04071 |  -0.0002956|   9.971503|  0.0000001|  0.9997529|  0.9999488|
| Niben101Scf04099g08003 |  -0.0000882|   8.674557|  0.0000000|  0.9999160|  0.9999727|
| Niben101Scf18308g01001 |   0.0000629|   6.922904|  0.0000000|  0.9999312|  0.9999727|
