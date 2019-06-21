##' compute_rank_score
##'
##'
##' @title psSubpathway internal functions
##' @description Compute rank score.
##' @usage SubSEA
##' @usage DCSA

##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong

compute_rank_score <- function(sort_idx_vec,n.genes){
  tmp <- rep(0, n.genes)
  tmp[sort_idx_vec] <- abs(seq(from=n.genes,to=1) - n.genes/2)
  return (tmp)
}

