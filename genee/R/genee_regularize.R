#' genee test using regularized effect size estimators.
#'
#' @param betas_ols A vector of effect size estimators (OLS betas).
#' @param alpha A tuning parameter determining the type of regularization. More specifically, alpha = 0: Ridge Regression; 0 < alpha < 1: Elastic Net; alpha = 1: LASSO.
#' @param ld A matrix of pairwise correlations statistics (i.e., r).
#' @param prior_weight A vector specifying the prior weight for each SNP.
#' @param gene_list A list where each element represents a vector containing the indices of all the SNPs in one set (e.g. a set could refer to a gene).
#' @return A list where the first element is a vector of test statistics and the second element is a vector of p values.
#' @export
#' @examples
#' x1 = c(0, 1, 0)
#' x2 = c(0, 1, 2)
#' x3 = c(1, 1, 2)
#' x = cbind(x1 ,x2, x3)
#' ld = cor(x)
#' betas = c(2, 3, -2)
#' prior_weight = c(1, 1, 1)
#' gene_list=list()
#' gene_list[[1]]=c(1, 2)
#' gene_list[[2]]=c(2, 3)
#' genee_regularize(betas, 0.5, ld, prior_weight, gene_list)
genee_regularize <- function(betas_ols, alpha, ld, prior_weight, gene_list){

  #run regularized regression
  beta_regularized = genee_regularize_regression(betas = betas_ols, ld = ld, alpha = alpha)

  #run EM to derive epsilon_effects
  epsilon_effect = genee_EM(betas = beta_regularized)

  #run test
  tempresults = genee_loop(betas = beta_regularized, ld = ld, epsilon_effect = epsilon_effect, prior_weight = prior_weight, gene_list = gene_list)

  #return results
  return(tempresults)
}
