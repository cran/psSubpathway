##' ssgsea
##'
##'
##' @title  psSubpathway internal functions
##' @description Single sample GSEA calculates a gene set enrichment score per sample as the normalized difference in
##' empirical cumulative distribution functions of gene expression ranks inside and outside the gene set.
##' @usage SubSEA
##' @usage DCSA

##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong
ssgsea <- function(X, geneSets, alpha=0.25, normalization=TRUE) {

  p <- nrow(X)
  n <- ncol(X)
  R <- apply(X, 2, function(x,p) as.integer(rank(x)), p)


  es <- sapply(1:n, function(j, R, geneSets, alpha) {

    geneRanking <- order(R[, j], decreasing=TRUE)
    es_sample <- NA
      es_sample <- sapply(geneSets, rndWalk, geneRanking, j, R, alpha)
    unlist(es_sample)
  }, R, geneSets, alpha)

  if (length(geneSets) == 1)
    es <- matrix(es, nrow=1)

  if (normalization) {

    es <- apply(es, 2, function(x, es) x / (range(es)[2] - range(es)[1]), es)
  }

  if (length(geneSets) == 1)
    es <- matrix(es, nrow=1)

  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)


  return(es)
}
