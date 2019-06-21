##' getgenedensity
##'
##'
##' @title psSubpathway internal functions
##' @description Calculating the Kernel estimation of gene.
##' @importFrom stats ecdf
##'
##' @usage SubSEA
##' @usage DCSA
##' @useDynLib psSubpathway

##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong

getgenedensity <- function(expr, sample.idxs, rnaseq=FALSE, kernel=TRUE){
  n.test.samples <- ncol(expr)
  n.genes <- nrow(expr)
  n.density.samples <- length(sample.idxs)

  gene.density <- NA
  if (kernel) {
    A = .C("matrix_density_R",
           as.double(t(expr[ ,sample.idxs, drop=FALSE])),
           as.double(t(expr)),
           R = double(n.test.samples * n.genes),
           n.density.samples,
           n.test.samples,
           n.genes,
           as.integer(rnaseq))$R

    gene.density <- t(matrix(A, n.test.samples, n.genes))
  } else {
    gene.density <- t(apply(expr, 1, function(x, sample.idxs) {
      f <- ecdf(x[sample.idxs])
      f(x)
    }, sample.idxs))
    gene.density <- log(gene.density / (1-gene.density))
  }

  return(gene.density)
}
