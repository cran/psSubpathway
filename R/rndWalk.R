##' rndWalk
##'
##'
##' @title psSubpathway internal functions
##' @description Calculating random walks.
##' @usage SubSEA
##' @usage DCSA

##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong

rndWalk <- function(gSetIdx, geneRanking, j, R, alpha) {
  indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
  stepCDFinGeneSet <- cumsum((abs(R[geneRanking, j]) *
                                indicatorFunInsideGeneSet)^alpha) /
    sum((abs(R[geneRanking, j]) *
           indicatorFunInsideGeneSet)^alpha)
  stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet) /
    sum(!indicatorFunInsideGeneSet)
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet

  sum(walkStat)
}
