##' ks_test_m
##'
##'
##' @title psSubpathway internal functions
##' @description Calculating subpathway Variation score.
##' @usage SubSEA
##' @usage DCSA
##' @useDynLib psSubpathway
##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong
##' @export

ks_test_m <- function(gset_idxs, gene.density, sort.idxs, mx.diff=TRUE,
                      abs.ranking=FALSE, tau=1){

  n.genes <- nrow(gene.density)
  n.samples <- ncol(gene.density)
  n.geneset <- length(gset_idxs)

  geneset.sample.es = .C("ks_matrix_R",
                         as.double(gene.density),
                         R = double(n.samples),
                         as.integer(sort.idxs),
                         n.genes,
                         as.integer(gset_idxs),
                         n.geneset,
                         as.double(tau),
                         n.samples,
                         as.integer(mx.diff),
                         as.integer(abs.ranking))$R

  return(geneset.sample.es)
}
