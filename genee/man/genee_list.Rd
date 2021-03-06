% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genee_list.R
\name{genee_list}
\alias{genee_list}
\title{Derive gene list using hg19 genes from UCSC browser.}
\usage{
genee_list(glist.hg19, all_chr, all_pos, upper, lower)
}
\arguments{
\item{glist.hg19}{A table containing all the information for hg19 genes from UCSC genome browser.}

\item{all_chr}{A vector specifying the chromosome for each SNP.}

\item{all_pos}{A vector specifying the physical location for each SNP.}

\item{upper}{A number specifying the number of base pairs downstream of the gene to be included.}

\item{lower}{A number specifying the number of base pairs upstream of the gene to be included.}
}
\value{
A list where the first element is a table containing information of all the genes included in the analysis, and the second element gives a list of indices for each gene.
}
\description{
Derive gene list using hg19 genes from UCSC browser.
}
