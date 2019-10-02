##' ssgsea
##'
##'
##' @title psSubpathway internal functions
##' @description Single sample GSEA (ssGSEA) calculates a gene set enrichment score
##' per sample as the normalized difference in empirical cumulative distribution
##' functions of gene expression ranks inside and outside the gene set.
##'
##' @usage SubSEA
##' @usage DCSA
##' @useDynLib psSubpathway

##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong

ssgsea <- function(X, geneSets, alpha=0.25, parallel.sz,
                   parallel.type, normalization=TRUE) {

  p <- nrow(X)
  n <- ncol(X)


  R <- apply(X, 2, function(x,p) as.integer(rank(x)), p)

  haveParallel <- isPackageLoaded("parallel")
  haveSnow <- isPackageLoaded("snow")

  cl <- makeCl <- parSapp <- stopCl <- mclapp <- detCor <- nCores <- NA
  if (parallel.sz > 1 || haveParallel) {

    if (!haveParallel) {  ## use snow
      makeCl <- get("makeCluster", mode="function")
      parSapp <- get("parSapply", mode="function")
      stopCl <- get("stopCluster", mode="function")


      cl <- makeCl(parallel.sz, type = parallel.type)
    } else {             ## use parallel

      mclapp <- get('mclapply', envir=getNamespace('parallel'))
      detCor <- get('detectCores', envir=getNamespace('parallel'))
      nCores <- detCor()
      options(mc.cores=nCores)
      if (parallel.sz > 0 && parallel.sz < nCores)
        options(mc.cores=parallel.sz)

    }
  }

  es <- sapply(1:n, function(j, R, geneSets, alpha) {
    geneRanking <- order(R[, j], decreasing=TRUE)
    es_sample <- NA
    if (parallel.sz == 1 || (is.na(cl) && !haveParallel))
      es_sample <- sapply(geneSets, rndWalk, geneRanking, j, R, alpha)
    else {
      if (is.na(cl))
        es_sample <- mclapp(geneSets, rndWalk, geneRanking, j, R, alpha)
      else
        es_sample <- parSapp(cl, geneSets, rndWalk, geneRanking, j, R, alpha)
    }

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

  es
}

