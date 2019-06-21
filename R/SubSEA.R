##' SubSEA
##'
##'
##' @title Subtype Set Enrichment Analysis (SubSEA).
##' @description The SubSEA (Subtype Set Enrichment Analysis) method to mine the specific subpathways of each sample
##' Subtype.
##' @param expr Matrix of gene expression values (rows are genes, columns are samples).
##' @param input.cls Input sample class vector (phenotype) file in CLS format.
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
##' This function uses the GSVA method to calculate the subpathway activity profile based on the gene expression
##' profile and subpathway list.Then we calculate the sample enrichment score (SES) of each subpathway by Subtype Set
##' Enrichment Analysis (SubSEA).We permute the gene labels and recompute the SES for the permuted data, which generates
##' a null distribution for the SES.The P value and the fdr value are calculated according to the perturbation analysis.
##'
##' @return A list containing the results of the SubSEA and the subpathway activity profile.
##' @author Xudong Han,
##' Junwei Han,
##' Qingfei Kong
##' @examples
##' # load depend package.
##' require(GSVA)
##' require(parallel)
##' # get breast cancer disease subtype gene expression profile.
##' Bregenematrix<-get("Subgenematrix")
##' # get path of the sample disease subtype files.
##' Subtypelabels<- system.file("extdata", "Sublabels.cls", package = "psSubpathway")
##' \donttest{SubSEA(Bregenematrix,input.cls=Subtypelabels,nperm=50,fdr.th=0.01,parallel.sz=2)}
##' # get the result of the SubSEA function
##' SubSEAresult<-get("Subspwresult")
##' str(SubSEAresult)
##' head(SubSEAresult$Basal)
##'
##' # Simulated gene matrix
##' genematrix <- matrix(rnorm(500*40), nrow=500, dimnames=list(1:500, 1:40))
##' # Construct subpathway list data.
##' subpathwaylist <- as.list(sample(2:100, size=20, replace=TRUE))
##' subpathwaylist <- lapply(subpathwaylist, function(n) sample(1:500, size=n, replace=FALSE))
##' names(subpathwaylist)<-c(paste(rep("spw",20),c(1:20)))
##' # Construct sample labels data.
##' subtypelabel<-list(phen=c("subtype1","subtype2","subtype3","subtype4"),
##'                    class.labes=c(rep("subtype1",10),rep("subtype2",10),
##'                    rep("subtype3",10),rep("subtype4",10)))
##' SubSEAcs<-SubSEA(genematrix,subtypelabel,subpathwaylist,nperm=10,parallel.sz=1)
##' str(SubSEAcs)
##'
##' @importFrom parallel parLapply
##' @importFrom parallel detectCores
##' @importFrom parallel makeCluster
##' @importFrom parallel clusterExport
##' @importFrom parallel clusterEvalQ
##' @importFrom parallel stopCluster
##' @importFrom GSVA gsva
##' @importFrom GSVA filterGeneSets
##' @importFrom stats p.adjust
##' @importFrom stats na.omit
##' @useDynLib psSubpathway,.registration = TRUE
##' @export



SubSEA<-function(expr,input.cls="",subpathwaylist="Symbol",kcdf="Gaussian",
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
  spwnames<-names(spwlist)

  if(method=="gsva"){
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
        cl <- makeCluster(num_workers)
        clusterExport(cl,"ks_test_m")


        clusterEvalQ(cl,library(GSVA))

        rdeslist<-parLapply(cl,1:nperm,function(n,zgindex,n.gset,spwlist,rank.scores,sort.sgn.idxs,mx.diff){

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
        stopCluster(cl)
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




  pn<-length(phen)

  if(parallel.sz==0){
    num_workers <- parallel::detectCores()
  }else{
    num_workers<-parallel.sz
  }
  cl2 <- makeCluster(num_workers)
  clusterExport(cl2,"FastSEAscore")
  rddata<-parLapply(cl2,1:pn,function(p,rdeslist,phen,sampleV,nperm){
    phen1<-phen[p]
    samlabel<-sampleV
    samlabel[samlabel!=phen1]<-0
    samlabel[samlabel==phen1]<-1
    ident<-as.numeric(samlabel)

    rdESmatrix<-lapply(1:nperm, function(x){
      rdm<-rdeslist[[x]]
      ryes<-apply(rdm,1,function(hang){
        pxindex<-order(hang,decreasing = T)
        pxhang<-hang[pxindex]
        pxlab<-ident[pxindex]
        res<-FastSEAscore(pxlab,pxhang)
        return(res)
      })
      return(ryes)
    })

    rdESmatrix1<-do.call(cbind,rdESmatrix)
    gc()
    return(rdESmatrix1)
  },rdeslist,phen,sampleV,nperm)
  stopCluster(cl2)

  ESdata<-sapply(1:pn, function(p){
    phen1<-phen[p]
    samlabel<-sampleV
    samlabel[samlabel!=phen1]<-0
    samlabel[samlabel==phen1]<-1
    ident<-as.numeric(samlabel)

    es<-apply(spwmatrix, 1, function(hang){
      pxindex<-order(hang,decreasing = T)
      pxhang<-hang[pxindex]
      pxlab<-ident[pxindex]
      yes<-FastSEAscore(pxlab,pxhang)
      return(yes)
    })
    return(es)
  })

  reportlist<-lapply(1:pn,function(j){
    ESvalue<-ESdata[,j]
    prdm<-rddata[[j]]

    p.vals <- NULL
    for(i in 1:length(ESvalue)){
      if(ESvalue[i]>=0){
        p.vals[i]<-sum(prdm[i,] >= ESvalue[i])/sum(prdm[i,]>=0)
      }else{
        p.vals[i]<-sum(prdm[i,] < ESvalue[i])/sum(prdm[i,]<0)

      }
    }
    ESscore<-ESvalue
    fdr<-p.adjust(p.vals,"BH",length(p.vals))
    sxindex<-which(fdr<=fdr.th)

    if(is.list(subpathwaylist)==TRUE){
      spwtitle<-names(subpathwaylist)
    }else{
      spwtitle<-get("spwtitle")
    }

    gene1<- genetag[sxindex]
    spwtitle1<-spwtitle[sxindex]

    result<-data.frame(SubpathwayID=names(spwlist)[sxindex],PathwayName=spwtitle1,SubpathwayGene=gene1,SES=ESscore[sxindex],Pvalue=p.vals[sxindex],FDR=fdr[sxindex])
    return(result)
  })
  reportlist[[length(phen)+1]]<-spwmatrix
  names(reportlist)<-c(phen,"spwmatrix")
  return(reportlist)

}
