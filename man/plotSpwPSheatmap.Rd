% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSpwPSheatmap.R
\name{plotSpwPSheatmap}
\alias{plotSpwPSheatmap}
\title{Plot subpathway phenotypic significant heat map}
\usage{
plotSpwPSheatmap(inputdata, spwid = "")
}
\arguments{
\item{inputdata}{A list of result data generated by function `SubSEA` or `DCSA`.}

\item{spwid}{The subpathway id which the user wants to plot.}
}
\value{
a plot
}
\description{
Visualize subpathway activity significant heat map between phenotypes
}
\details{
plotSpwPSheatmap


A heat map of the significance P-value of the activity of the subpathway between the phenotypes.
The rows and columns of the heat map are sample phenotype labels. The values shown in the heat map are the T-test
P values between the activity of the subpathway corresponding to the two phenotypes. The lower the number in the cells
in the heat map, the greater the difference in the activity of the subpathways between the two phenotypes.
}
\examples{
# get the Subspwresult which is the result of SubSEA method.
Subspwresult<-get("Subspwresult")
# plot significant heat map between the activity of the subpathway in each subtype of breast cancer.
plotSpwPSheatmap(Subspwresult,spwid="00120_9")
# get the DCspwresult which is the result of DCSA method.
DCspwresult<-get("DCspwresult")
# plot significant heat map between the activity of the subpathway in each stage of breast cancer.
plotSpwPSheatmap(DCspwresult,spwid="00982_2")
}
\author{
Xudong Han,
Junwei Han,
Qingfei Kong
}
