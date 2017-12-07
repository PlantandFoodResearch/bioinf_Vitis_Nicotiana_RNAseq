
## Load edgeR
if(!require(edgeR)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  library(edgeR)
}

filename <- '/workspace/hradxj/karmun_awesome_experiment/010.edgeR_Nb/GRLaV3_Nb_EdgeR.tab'

GRLaV3_vs_Nb <- read.delim(filename, header=TRUE)
dim(GRLaV3_vs_Nb)
#GRLaV3_vs_Nb2 <- read.table(filename, header=TRUE, sep="\t")

## Sanity check names problem
## Missing Gene name - Looks like an off by one error indroduced in preprocessing/filtering
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% ""),]

## Duplicated names
badnames <- names(table(GRLaV3_vs_Nb$Gene))[table(GRLaV3_vs_Nb$Gene)>1]
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% badnames),]

# Row names are not unique for some reason. Crap annotations again.
# Assign unique row names using the preexisting names, but making them
# unique where required by adding ".1, .2" etc.
rownames(GRLaV3_vs_Nb) <- make.names(GRLaV3_vs_Nb[,1], unique=TRUE)

## Sanity check
GRLaV3_vs_Nb[which(GRLaV3_vs_Nb$Gene %in% badnames),]

# Remove the redundant column where the row names came from
GRLaV3_vs_Nb[,1] <- NULL


panel.points <- function(x, y, col=NULL, loess = FALSE,  maTransform = TRUE, ...){
  if(is.null(col)) {
    dcols <- densCols(x,y)
    if(maTransform) {
      points((x+y)/2, y-x, col=dcols, ...)
    } else {
      points(x, y, col=dcols, ...)
    }
  } else {
    if(maTransform) {
      points((x+y)/2, y-x, col=col, ...)
    } else {
      points(x, y, col=col, ...)
    }
  }
  if(loess) {
    loessObj <- loess(y~x)
    yhat <- fitted(loessObj)
    o <- order(yhat)
    lines(x[o], yhat[o], col="red")
  }
  if(maTransform) {
    abline(h=0, cex=0.01, col="darkgrey")
  } else {
    abline(0,1, cex=0.01, col="darkgrey")
  }
}


panel.corr <- function(x,y, maTransform = TRUE, ...){
  NAs <- (is.na(x) | is.na(y))
  correlation <- round(cor(x[!NAs], y[!NAs]), 2)
  if(maTransform) {
    text(x=10, y=0, labels=substitute(rho==correlation, list(correlation=correlation)))
  } else {
    text(x=7.5, y=7.5, labels=substitute(rho==correlation, list(correlation=correlation)))
  }
}

## Healthy.2 replicate is different from the other samples (including the other biological replicate)
pairs(log2(GRLaV3_vs_Nb+0.5), labels=colnames(tmp), upper.panel=panel.points, lower.panel=panel.corr, 
      pch=16, cex=0.4, loess=FALSE, maTransform=FALSE)

pairs(log2(GRLaV3_vs_Nb+0.5), labels=colnames(tmp), upper.panel=panel.points, lower.panel=panel.corr, 
      xlim=c(0,20), ylim=c(-8,8), pch=16, cex=0.4, loess=FALSE)

plotMD.DGEList(de_genes)

plotMDS(log2(GRLaV3_vs_Nb+0.5), top=1000)

# set the replicate relationship.. columns 1 and 2 are reps of each other
## Healthy vs Infected
replicate_relationship <- factor(c(1,1,2))
# construct the DGEList object
de_genes <- DGEList(counts=GRLaV3_vs_Nb,group=replicate_relationship)

# Calculate normalisation factors
de_genes <- calcNormFactors(de_genes)
# Examine normalisation factors as a sanity check
de_genes$samples
# We can see that Benth.Healthy2 had a lot less sequencing.

# Nonsensical results can happen if we do not remove genes
# that are not biologically relevant because of low expression
# Filter out low expressed genes on the basis of the count per million
# This filters out genes that are less than cpmlimit, and 
# that are not expressed in libexplimit number of libraries
cpmlimit    <- 1000
libexplimit <- 3
keep <- rowSums(cpm(de_genes)>cpmlimit) >= libexplimit


tmp <- log2(cpm(de_genes))
par(pty="s")
pairs(tmp, pch=16, cex=0.3, asp=1)


## Sanity check how many are kept
table(keep)


de_genes_kept <- de_genes[keep, , keep.lib.sizes=FALSE]
# Calculate normalisation factors again
de_genes_kept <- calcNormFactors(de_genes)
# Examine normalisation factors as a sanity check again
de_genes_kept$samples


# set experiment design
design <- model.matrix(~replicate_relationship)
# estimate dispersion
de_genes_kept <- estimateDisp(de_genes_kept, design)




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
fitglm <- glmFit(de_genes_kept, design)
lrtglm <- glmLRT(fitglm, coef=2)
topTags(lrtglm, n=10)


# Use the quasi-likelihood model to obtain DE genes.
# From edgeR manual: 
#"While the likelihood ratio test is a more obvious choice for inferences with GLMs, the QL
#F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It
#provides more robust and reliable error rate control when the number of replicates is small"
fitqlm <- glmQLFit(de_genes_kept, design)
qlftest <- glmQLFTest(fitqlm, coef=2)
topTags(qlftest, n=10)
