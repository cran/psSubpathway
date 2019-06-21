##' DCSA
##'
##'
##' @title Dynamic Changed Subpathway Analysis (DCSA)
##' @description This function will perform the Dynamic Changed Subpathway Analysis (DCSA) method to estimate the dynamic
##' changed subpathways associated with the sample phenotypes (like developmental stage of cancer).
##' @param expr Matrix of gene expression values (rows are genes, columns are samples).
##' @param input.cls Input sample phenotype class vector file in CLS format.
##' @param subpathwaylist Character string denoting the gene label of the subpahtway list is `Entrezid` (default) or `Symbol`.
##' Users can also enter their own subpathway list data. This list should be consistent with the gene label in the input gene expression profile.
##' @param kcdf Character string denoting the kernel to use during the non-parametric estimation of the cumulative
##' distribution function of expression levels across samples when method="gsva". By default, `kcdf="Gaussian"` which
##' is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale,
##' RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as those derived from
##' RNA-seq experiments, then this argument should be set to `kcdf="Poisson"`.
##' @param method Method to employ in the estimation of subpathway enrichment scores per sample.By default this is set to
##' `gsva` (Hänzelmann et al, 2013) and other options are `ssgsea` (Barbie et al, 2009).
##' @param min.sz Minimum size of the resulting subpathway.
##' @param max.sz Maximum size of the resulting subpathway.
##' @param nperm Number of random permutations (default: 100).
##' @param fdr.th Cutoff value for fdr. Only subpathway with lower fdr.th are listed (default: 1).
##' @param mx.diff Offers two approaches to calculate the sample enrichment score (SES) from the KS random walk statistic.
##' `mx.diff=FALSE`: SES is calculated as the maximum distance of the random walk from 0. `mx.diff=TRUE` (default): SES is
##' calculated as the magnitude difference between the largest positive and negative random walk deviations.
##' @param parallel.sz Number of processors to use when doing the calculations in parallel. If this argument is left with its
##' default value (parallel.sz=0) then it will use all available core processors unless we set this argument with a smaller number.
##' @details
##' This function constructs a subpathway activity profile  by using subpathway list and gene expression
##' profiling as input to the GSVA method. Next,we used the information theoretic measure of statistical dependence,
##' mutual information (MI), to estimate the dynamic changed subpathways associated with the sample phenotypes.
##' Finally we used the perturbation analysis of the gene label rearrangement to estimating the statistical significance.
##'
##' @return A list containing the results of DCSA and subpathway activity profile.
##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong
##' @examples
##' # load depend package.
##' require(GSVA)
##' require(parallel)
##' require(mpmi)
##' # get ACC disease stage gene expression profiling.
##' ACCgenematrix<-get("DCgenematrix")
##' # get path of the sample disease stage phenotype files.
##' Stagelabels<-system.file("extdata", "DClabels.cls", package = "psSubpathway")
##' # perform the DCSA method.
##' \donttest{DCSA(ACCgenematrix,input.cls=Stagelabels,nperm=50,fdr.th=0.01,parallel.sz=2)}
##' # get the result of the SubSEA function
##' DCSAresult<-get("DCspwresult")
##' str(DCSAresult)
##' head(DCSAresult$DCSA)
##'
##' # Simulated gene matrix.
##' genematrix <- matrix(rnorm(500*40), nrow=500, dimnames=list(1:500, 1:40))
##' # Construct subpathway list data.
##' subpathwaylist <- as.list(sample(2:100, size=20, replace=TRUE))
##' subpathwaylist <- lapply(subpathwaylist, function(n) sample(1:500, size=n, replace=FALSE))
##' names(subpathwaylist)<-c(paste(rep("spw",20),c(1:20)))
##' # Construct sample labels data.
##' stagelabel<-list(phen=c("stage1","stage2","stage3","stage4"),
##'                    class.labes=c(rep("stage1",10),rep("stage2",10),
##'                    rep("stage3",10),rep("stage4",10)))
##' DCSAcs<-DCSA(genematrix,stagelabel,subpathwaylist,nperm=10,parallel.sz=1)
##' str(DCSAcs)
##'
##' @importFrom parallel parLapply
##' @importFrom parallel detectCores
##' @importFrom parallel makeCluster
##' @importFrom parallel clusterExport
##' @importFrom parallel clusterEvalQ
##' @importFrom parallel stopCluster
##' @importFrom GSVA gsva
##' @importFrom GSVA filterGeneSets
##' @importFrom mpmi mminjk.pw
##' @importFrom stats na.omit
##' @importFrom stats p.adjust
##'
##' @useDynLib psSubpathway,.registration = TRUE
##' @export

DCSA<-function(expr,input.cls="",subpathwaylist="Symbol",kcdf="Gaussian",
               method="gsva",min.sz=1,max.sz=Inf,nperm=100,fdr.th=1,
               mx.diff=TRUE,parallel.sz=0){

  if(is.list(subpathwaylist)==TRUE){
    spwlist<-subpathwaylist
    genetag<-unlist(lapply(subpathwaylist,
                           function(x){
                             x<-as.character(x)
                             a<-paste(unlist(x),collapse=" ")
                             return(a)
                           }))
  }else{
    if(subpathwaylist=="Entrezid"){
      spwlist<-get("spwentrezidlist")
      genetag<-get("geneentrezid")
    }else if(subpathwaylist=="Symbol"){
      spwlist<-get("spwsymbollist")
      genetag<-get("genesymbol")
    }
  }

  if (kcdf == "Gaussian") {
    rnaseq <- FALSE
    kernel <- TRUE
  } else if (kcdf == "Poisson") {
    rnaseq <- TRUE
    kernel <- TRUE
  }


  allgene<-unlist(spwlist)
  qcindex<-match(rownames(expr),allgene)
  naindex<-which(is.na(qcindex)==TRUE)
  if(length(naindex)>0){
    expr<-expr[-naindex,]
  }else{
    expr<-expr
  }

  if(is.list(input.cls)) {
    CLS <- input.cls
  } else {
    CLS <- ReadClsFile(file=input.cls)
  }
  phen <- rev(CLS$phen)
  sampleV <-CLS$class.labes

  spwlist <- lapply(spwlist,
                    function(x ,y) na.omit(match(x, y)),
                    rownames(expr))
  spwlist <- filterGeneSets(spwlist,
                            min.sz=max(1, min.sz),
                            max.sz=max.sz)


  if(method=="gsva"){

    spwnames<-names(spwlist)
    n.samples <- ncol(expr)
    n.genes <- nrow(expr)
    zgindex<-seq(1,n.genes)
    sample.idxs<-c(1:n.samples)

    gene.density <- getgenedensity(expr, sample.idxs, rnaseq, kernel)

    rank.scores <- rep(0, n.genes)
    sort.sgn.idxs <- apply(gene.density, 2, order, decreasing=TRUE)
    rank.scores <- apply(sort.sgn.idxs, 2, compute_rank_score,n.genes)

    n.gset<-length(spwlist)

    spwmatrix<-t(sapply(spwlist, ks_test_m, rank.scores, sort.sgn.idxs,
                        mx.diff=mx.diff, abs.ranking=FALSE,
                        tau=1))
    row.names(spwmatrix)<-spwnames
    colnames(spwmatrix)<-sampleV

    if(parallel.sz==1){
      rdeslist<-lapply(1:nperm, function(n,zgindex,n.gset,spwlist,rank.scores,sort.sgn.idxs,mx.diff){

        rdgs.list<-lapply(1:n.gset, function(s){
          ng<-length(spwlist[[s]])
          rdgs<-sample(zgindex,ng,replace = F)
          return(rdgs)
        })

        m <- t(sapply(rdgs.list, ks_test_m, rank.scores, sort.sgn.idxs,
                      mx.diff=mx.diff, abs.ranking=FALSE,
                      tau=1))

        return(m)

      },zgindex,n.gset,spwlist,rank.scores,sort.sgn.idxs,mx.diff)

    }else{
      if(parallel.sz==0){
        num_workers <- parallel::detectCores()
      }else{
        num_workers<-parallel.sz
      }
      cl3 <- makeCluster(num_workers)
      clusterExport(cl3,"ks_test_m")

      clusterEvalQ(cl3,library(GSVA))

      rdeslist<-parLapply(cl3,1:nperm,function(n,zgindex,n.gset,spwlist,rank.scores,sort.sgn.idxs,mx.diff){

        rdgs.list<-lapply(1:n.gset, function(s){
          ng<-length(spwlist[[s]])
          rdgs<-sample(zgindex,ng,replace = F)
          return(rdgs)
        })

        m <- t(sapply(rdgs.list, ks_test_m, rank.scores, sort.sgn.idxs,
                      mx.diff=mx.diff, abs.ranking=FALSE,
                      tau=1))

        return(m)

      },zgindex,n.gset,spwlist,rank.scores,sort.sgn.idxs,mx.diff)
      stopCluster(cl3)
    }

  }
  if(method=="ssgsea"){

    genename<-row.names(expr)
    spwmatrix<-ssgsea(expr,spwlist,alpha=1)

    rdeslist<-lapply(1:nperm,function(n){
      expr1<-expr
      row.names(expr1)<-sample(genename,replace = F)

      m<-ssgsea(expr1,spwlist,alpha=1)
      return(m)
    })
  }


  hxx<-apply(spwmatrix,1,mminjk.pw,sampleV)

  if(parallel.sz==0){
    num_workers <- parallel::detectCores()
  }else{
    num_workers<-parallel.sz
  }
  cl4 <- makeCluster(num_workers)
  clusterEvalQ(cl4,library(mpmi))

  if(is.list(subpathwaylist)==TRUE){
    spwtitle<-names(subpathwaylist)
  }else{
    spwtitle<-get("spwtitle")
  }
  rdh<-parLapply(cl4,1:nperm,function(n,rdeslist,sampleV){
    rdspwmatrix<-rdeslist[[n]]
    yrdhxx<-apply(rdspwmatrix,1,function(hang){
      hxx1<-mminjk.pw(hang,sampleV)
      return(hxx1)
    })
    return(yrdhxx)
  },rdeslist,sampleV)
  stopCluster(cl4)
  rdhxx<-do.call(cbind,rdh)


  pz<-NULL
  for(i in 1:length(rdhxx[,1])){

    pz[i]<-sum(rdhxx[i,] >= hxx[i])/nperm

  }
  fdr<-p.adjust(pz,method = "BH",length(pz))
  sxindex<-which(fdr<=fdr.th)



  gene1<- genetag[sxindex]
  spwtitle1<-spwtitle[sxindex]

  report<-data.frame(SubpahwayID=names(spwlist)[sxindex],PahwayName=spwtitle1,SubpathwayGene=gene1,"D(M)"=hxx[sxindex],Pvalue=pz[sxindex],FDR=fdr[sxindex],
                     check.names = F)
  hxxreport<-list(report,spwmatrix)
  names(hxxreport)<-c("DCSA","spwmatrix")
  return(hxxreport)
}
