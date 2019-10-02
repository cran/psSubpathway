## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
require(psSubpathway)

## ------------------------------------------------------------------------
# load depend package.
 require(GSVA)
 require(parallel)
# get breast cancer disease subtype gene expression profile.
 Bregenematrix<-get("Subgenematrix")
# get path of the sample disease subtype files.
 Subtypelabels<- system.file("extdata", "Sublabels.cls", package = "psSubpathway")
# perform the SubSEA method.
#SubSEAresult<-SubSEA(Bregenematrix,
#                     input.cls=Subtypelabels,
#                     nperm=50,                     # Number of random permutations
#                     fdr.th=0.01,                  # Cutoff value for fdr.when fdr.th=1,get all subpathway
#                     parallel.sz=2)                # Number of processors to use  
 
# get the result of the SubSEA function
 SubSEAresult<-get("Subspwresult")
 str(SubSEAresult)
 head(SubSEAresult$Basal)

## ----message=FALSE, warning=FALSE----------------------------------------
# load depend package.
 require(mpmi)
# get ACC disease stage gene expression profiling.
 ACCgenematrix<-get("DCgenematrix")
# get path of the sample disease stage phenotype files.
 Stagelabels<-system.file("extdata", "DClabels.cls", package = "psSubpathway")
# perform the DCSA method.
#DCSAresult<-DCSA(ACCgenematrix,input.cls=Stagelabels,nperm=50,fdr.th=0.01,parallel.sz=2)
# get the result of the SubSEA function
 DCSAresult<-get("DCspwresult")
 str(DCSAresult)
 head(DCSAresult$DCSA)

## ----fig.height=3, fig.width=12------------------------------------------
# plot enrichment score curve of the subpathway 00120_9 in all breast cancer subtypes
 plotSubSEScurve(Subspwresult,
                 spwid="00120_9", # The subpathway id which subpahtway is ploted
                 phenotype="all")  # Which phenotypic phenotype set enrichment curve is drawn
                

## ----fig.height=4, fig.width=5-------------------------------------------

# plot enrichment score curve of the subpathway 00120_9 in the Basal breast cancer subtype.
 plotSubSEScurve(Subspwresult,spwid="00120_9",phenotype="Basal")

## ----fig.height=6, fig.width=10------------------------------------------
# plot the subpathway 00120_9 in the SubSEA function result.
 plotSpwACmap(Subspwresult,spwid="00120_9")
 

## ----fig.height=6, fig.width=10------------------------------------------
# plot the subpathway 00982_2 in the DCSA function result.
 plotSpwACmap(DCspwresult,spwid="00982_2")

## ----fig.height=8, fig.width=10------------------------------------------
# load depend package.
 require(pheatmap)
# plot significant up-regulation subpathway heat map specific for each breast cancer subtype.
 plotheatmap(Subspwresult,
             fdr.th=0.01,      # Cutoff value for fdr
             plotSubSEA=TRUE,  # Indicate that the input data is the result of the SubSEA function
             SES="positive",   # Obtain a subpathway with a positive SES value
             phenotype="all")


## ----fig.height=8, fig.width=10------------------------------------------
# plot significant down-regulation subpathway heat map specific for each breast cancer subtype.
 plotheatmap(Subspwresult,plotSubSEA=TRUE,fdr.th=0.01,SES="negative",phenotype="all")


## ----fig.height=8, fig.width=10------------------------------------------
# plot Basal subtype specific significant subpathway heat map.
 plotheatmap(Subspwresult,plotSubSEA=TRUE,fdr.th=0.01,SES="all",phenotype="Basal")


## ----fig.height=8, fig.width=10------------------------------------------
# plot a heat map of the subpathway that is significantly associated with breast cancer stages.
 plotheatmap(DCspwresult,plotSubSEA=FALSE,fdr.th=0.01)


## ----fig.height=6, fig.width=7-------------------------------------------
# get the Subspwresult which is the result of SubSEA method.
 Subspwresult<-get("Subspwresult")
 # plot significant heat map between the activity of the subpathway in each subtype of breast cancer.
 plotSpwPSheatmap(Subspwresult,spwid="00120_9")

## ----fig.height=6, fig.width=7-------------------------------------------
# get the DCspwresult which is the result of DCSA method. 
 DCspwresult<-get("DCspwresult")
##' # plot significant heat map between the activity of the subpathway in each stage of breast cancer.
 plotSpwPSheatmap(DCspwresult,spwid="00982_2")

## ----fig.height=6, fig.width=8, message=FALSE, warning=FALSE-------------
# load depend package.
 library(igraph)
# plot network graph of the subpathway 00982_2
plotSpwNetmap("00982_2")


