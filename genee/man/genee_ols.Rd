% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genee_ols.R
\name{genee_ols}
\alias{genee_ols}
\title{genee test using effect size estimators derived through ordinary least square method (no shrinkage).}
\usage{
genee_ols(betas_ols, ld, prior_weight, gene_list)
}
\arguments{
\item{betas_ols}{A vector of effect size estimators (OLS betas).}

\item{ld}{A matrix of pairwise correlations statistics (i.e., r).}

\item{prior_weight}{A vector specifying the prior weight for each SNP.}

\item{gene_list}{A list where each element represents a vector containing the indices of all the SNPs in one set (e.g., a set could refer to a gene).}
}
\value{
A list where the first element is a vector of test statistics; the second element is a vector of variances and the third element is a vector of p-values.
}
\description{
genee test using effect size estimators derived through ordinary least square method (no shrinkage).
}
\examples{
x1 = c(0, 1, 1)
x2 = c(0, 1, 2)
x3 = c(1, 1, 2)
x = cbind(x1 ,x2, x3)
ld = cor(x)
betas = c(0.1, 0.2, 0.1)
prior_weight = c(1, 1, 1)
gene_list=list()
gene_list[[1]]=c(1, 2)
gene_list[[2]]=c(2, 3)
genee_ols(betas, ld, prior_weight, gene_list)
}
