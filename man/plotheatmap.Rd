% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotheatmap.R
\name{plotheatmap}
\alias{plotheatmap}
\title{Plot a heatmap}
\usage{
plotheatmap(
  inputdata,
  plotSubSEA = TRUE,
  fdr.th = 1,
  SES = "positive",
  phenotype = "all"
)
}
\arguments{
\item{inputdata}{A list of result data generated by function `SubSEA` or `DCSA`.}

\item{plotSubSEA}{Determine the inputdata is the result data of function `SubSEA` (default:plotSubSEA=TRUE) or
`DCSA` (plotSubSEA=FLASE).}

\item{fdr.th}{Cutoff value for FDR. The only subpathway with lower FDR is plotted. (default: 1).}

\item{SES}{Parameter `SES` is useful only when `plotSubSEA` is TRUE. When `plotSubSEA=TRUE`,if `SES` is positive,
the subpathway with high-expression will be plotted. when it is negative, plot low-expression subpathways.}

\item{phenotype}{Parameter `phenotype` is useful only when `plotSubSEA` is TRUE. `phenotype` decides which phenotypic
significant subpathways to screen (which phenotypic result is applied to parameter `fdr.th` and `SES`.) and plot a
heat map of these subpathways.By default,`phenotype="all"` which will screen the subpathways of all phenotypes and plot
a heat map. When the user wants to plot a subpathway heat map of the specified phenotype, this parameter should be set
to the name of the phenotype.}
}
\value{
a heatmap
}
\description{
Plot a heatmap of subpathway activity profile based on the parameters set by the user.
}
\details{
plotheatmap


The subpathways are screened according to the conditions set by the user and a heat map of the activity of these subpathways is drawn.
}
\examples{
# load depend package.
library(pheatmap)
# get the Subspwresult which is the result of SubSEA function.
Subspwresult<-get("Subspwresult")
# get the DCspwresult which is the result of DCSA function.
DCspwresult<-get("DCspwresult")
# plot significant up-regulation subpathway heat map specific for each breast cancer subtype.
plotheatmap(Subspwresult,plotSubSEA=TRUE,fdr.th=0.01,SES="positive",phenotype="all")
# plot significant down-regulation subpathway heat map specific for each breast cancer subtype.
plotheatmap(Subspwresult,plotSubSEA=TRUE,fdr.th=0.01,SES="negative",phenotype="all")
# plot basal subtype specific significant subpathway heat map.
plotheatmap(Subspwresult,plotSubSEA=TRUE,fdr.th=0.01,SES="all",phenotype="Basal")
# plot adrenocortical cancer disease stages specific significant subpathway heat map.
plotheatmap(DCspwresult,plotSubSEA=FALSE,fdr.th=0.01)
}
\author{
Xudong Han,
Junwei Han,
Qingfei Kong
}
