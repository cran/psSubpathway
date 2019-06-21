##' SEAscore
##'
##'
##' @title psSubpathway internal functions
##' @description Get subtype set enrichment score and sample locations, etc.
##' @usage plotSubSEScurve
##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong


SEAscore<-function(labels.list,correl.vector = NULL){
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
  if (max.ES > - min.ES) {
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}
