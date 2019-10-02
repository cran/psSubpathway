\name{envData}
\alias{envData}
\alias{DCgenematrix}
\alias{DCspwresult}
\alias{spwtitle}
\alias{spwentrezidlist}
\alias{spwsymbollist}
\alias{spwnetworkdata}
\alias{Subgenematrix}
\alias{Subspwresult}
\alias{DClabes.cls}
\alias{Sublabels.cls}
\docType{data}

\title{The variables in the environment include subpathway list data, information of subpathway,an expression profile and a example result}
\description{
    We used the k-clique algorithm to divide the human pathway of the KEGG database into subpathways and eliminated the smaller module that had a overlap above 80%.
  These subpathway data will be divided into two parts: subpathway information(spwnetworkdata,spwtitle) and subpathway list data(spwentrezidlist: the gene Entrezid contained in each subpathway and spwsymbollist:the gene symbol contained in each subpathway).
  Subgenematrix is the gene expression profile of the breast cancer we selected from GDCTCGA and Subspwresult is the result of applying Subgenematrix to Subtype Set Enrichment Analysis (SubSEA).In order to reduce the memory, we delete the gene in the gene expression profile that is not in the gene of the subpathway list.
  We also selected the gene expression profile of the adrenocortical cancer (ACC) from GDCTCGA and performed Dynamic Changed Subpathway Analysis (DCSA). The data are DCgenematrix and DCspwresult
  In order to reduce the memory, we delete the gene in the gene expression profile that is not in the gene of the subpathway list.
  respectively. DClabes.cls and Sublabels.cls are the label vectors of the phenotypes of the samples of the two gene expression profiles.
  }

\format{An environment variable}
\author{Junwei Han\email{hanjunwei1981@163.com},Xudong Han \email{HanXD1217@163.com}}

\keyword{data}
