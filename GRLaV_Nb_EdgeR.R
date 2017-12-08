
library(edgeR)
library(limma)

rm(list=ls())


## Graphics functions for pairs plots
source("functions.R")

filename <- '/workspace/hradxj/karmun_awesome_experiment/010.edgeR_Nb/GRLaV3_Nb_EdgeR.tab'

GRLaV3_vs_Nb <- read.delim(filename, header=TRUE)
dim(GRLaV3_vs_Nb)
#GRLaV3_vs_Nb2 <- read.table(filename, header=TRUE, sep="\t")

## Sanity check names problem
## Missing Gene name - Looks like an off by one error indroduced in preprocessing/filtering
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% ""), ]

## Duplicated names
badnames <- names(table(GRLaV3_vs_Nb$Gene))[table(GRLaV3_vs_Nb$Gene)>1]
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% badnames), ]

# Row names are not unique for some reason. Crap annotations again.
# Assign unique row names using the preexisting names, but making them
# unique where required by adding ".1, .2" etc.
rownames(GRLaV3_vs_Nb) <- make.names(GRLaV3_vs_Nb[,1], unique=TRUE)

## Sanity check
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% badnames),]

# Remove the redundant column where the row names came from
GRLaV3_vs_Nb[,1] <- NULL

## Healthy.2 replicate is different from the other samples (including the other biological replicate)

par(pty="s")

pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr,
      pch=16, cex=0.4, asp=1)

pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr,
      maTransform=TRUE, xlim=c(0,20), ylim=c(-8,8), pch=16, cex=0.4, asp=1)


plotMDS(log2(GRLaV3_vs_Nb+0.5), top=1000)


## How many samples are 0 in All samples?
# Filtering these has no effect on cpm/rpkm calculations
apply(GRLaV3_vs_Nb==0, 1, all)  %>%  table


# set the replicate relationship.. columns 1 and 2 are reps of each other
## Healthy vs Infected
replicate_relationship <- factor(c(1,1,2))
# construct the DGEList object
all_genes <- DGEList(counts=GRLaV3_vs_Nb, group=replicate_relationship)

for( i in seq(ncol(all_genes))) {
  plotMD.DGEList(all_genes, column=i)
  abline(h=0, col="grey", lty=2)
}


breaks <- seq(0, ceiling(max(log2(all_genes$counts))))
tweak  <- 1

## Distribution of counts
for(i in colnames(all_genes)) {
  hist(log2(all_genes$counts+tweak)[, i], breaks=breaks, main=i, xlab="log2(counts)")
}

# Calculate normalisation factors
all_genes <- calcNormFactors(all_genes, method="TMM")

# Examine normalisation factors as a sanity check
all_genes$samples
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

## Sanity check how many are kept
table(keep)


## Visualize filtering
pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr,
      pch=16, cex=0.1, asp=1, col=c("black","grey")[(keep)+1], main="Filtering effect")

pairs(log2(GRLaV3_vs_Nb), labels=colnames(GRLaV3_vs_Nb), upper.panel=panel.points, lower.panel=panel.corr,
      maTransform=TRUE, xlim=c(0,20), ylim=c(-8,8), pch=16, cex=0.1, asp=1, col=c("black", "grey")[(keep)+1], main="MAplot Filtering effect")


## Distribution of counts after filtering
for(i in colnames(all_genes)) {
  hist(log2(all_genes$counts+tweak)[keep , i], breaks=breaks, main=i, xlab="log2(counts)")
}


## Alternative filter which overrides line:89
keep <- !apply(GRLaV3_vs_Nb==0, 1, all)

## Sanity check how many are kept
table(keep)

## Must use which here
all_genes_kept <- all_genes[which(keep), keep.lib.sizes=FALSE]

## Sanity check
dim(all_genes)
dim(all_genes_kept)

# Calculate normalisation factors again
all_genes_kept <- calcNormFactors(all_genes_kept, method="TMM")

## Sanity check
apply(GRLaV3_vs_Nb,2, sum)
calcNormFactors(all_genes, method="TMM")
all_genes_kept$samples

## Benth.Healthy.2 is an order of magnitude smaller!
apply(GRLaV3_vs_Nb,2, sum)[1] / apply(GRLaV3_vs_Nb,2, sum)[2]

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


threshold <- 0.05
de_genes  <- topN$table[topN$table$FDR<threshold,]

## DE genes
nrow(de_genes)

## First 10
head(de_genes, n=10)

## estimate pi0
pvals <- topN$table$PValue

hist(pvals)

## Investigating abnormal right hand peak in pvalue histogram
# ind <- which(pvals>=1)
# rownames(topN$table)[ind]
# GRLaV3_vs_Nb[rownames(topN$table)[ind],]

pi0_hat <- limma::convest(pvals)
print(pi0_hat)

# Use the quasi-likelihood model to obtain DE genes.
# From edgeR manual: 
#"While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small"
fitqlm  <- glmQLFit(all_genes_kept, design)
qlftest <- glmQLFTest(fitqlm, coef=2)
topN2   <- topTags(qlftest, n=nrow(lrtglm))

threshold <- 0.05
de_genes2 <- topN2$table[topN2$table$FDR<threshold,]

## DE genes
nrow(de_genes2)
## First 10
head(de_genes2, n=10)


## estimate pi0
pvals <- topN$table$PValue

hist(pvals)

## Investigating abnormal right hand peak in pvalue histogram
# ind <- which(pvals>=1)
# rownames(topN$table)[ind]
# GRLaV3_vs_Nb[rownames(topN$table)[ind],]

pi0_hat <- limma::convest(pvals)
print(pi0_hat)
