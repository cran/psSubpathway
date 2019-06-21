##' FastSEAscore
##'
##'
##' @title psSubpathway internal functions
##' @description Fast calculate phenotypic set enrichment score.
##' @usage SubSEA
##' @usage DCSA
##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong
##' @export

FastSEAscore<-function(labels.list,correl.vector = NULL){
  tag.indicator <- labels.list
  no.tag.indicator <- 1 - tag.indicator
  N <- length(labels.list)
  Nh <- length(labels.list[labels.list==1])
  Nm <-  N - Nh
  correl.vector <- abs(correl.vector)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)

  max.ES <- max(RES)
  min.ES <- min(RES)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  return(ES)
}
