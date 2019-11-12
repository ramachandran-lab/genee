#' genee test using regularized effect size estimators.
#'
#' @param betas_ols A vector of effect size estimators (OLS betas).
#' @param alpha A tuning parameter determining the type of regularization. More specifically, alpha = 0: Ridge Regression; 0 < alpha < 1: Elastic Net; alpha = 1: LASSO.
#' @param ld A matrix of pairwise correlations statistics (i.e., r).
#' @param ld_path Path to the ld matrix memory mapping file (.desc, .bin) including the prefix of the file name. Make sure that .desc and .bin are in the same directory.
#' @param prior_weight A vector specifying the prior weight for each SNP.
#' @param gene_list A list where each element represents a vector containing the indices of all the SNPs in one set (e.g. a set could refer to a gene).
#' @param nfolds A number indicates how many folds to do for cross validation for regularize regression. The default is 10. If it's fed with 0, then it will perform no cross validation and choose lambda with smallest penalty.
#' @param lambda.min The smallest lambda to use for regularize regression, as a fraction of the largest lambda. Default is 0.0001.
#' @param nlambda The number of lambda values to use. Default is 100.
#' @return A list where the first element is a vector of test statistics and the second element is a vector of p values.
#' @export
#' @examples
#' x1 = c(0, 1, 0)
#' x2 = c(0, 1, 2)
#' x3 = c(1, 1, 2)
#' x = cbind(x1 ,x2, x3)
#' ld = cor(x)
#' write.table(ld, file = "myld", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
#' ld_name = "myld"
#' prepare_LD(ld_name)
#' ld_path = "myld"
#' betas = c(2, 3, -2)
#' prior_weight = c(1, 1, 1)
#' gene_list=list()
#' gene_list[[1]]=c(1, 2)
#' gene_list[[2]]=c(2, 3)
#' genee_regularize_bigdata(betas, 0.5, ld_path, prior_weight, gene_list, 5, 0.001, 50)
genee_regularize_bigdata <- function(betas_ols, alpha, ld, ld_path, prior_weight, gene_list, nfolds, lambda.min, nlambda){

  #run regularized regression
  beta_regularized = genee_regularize_regression_bigdata(betas = betas_ols, ld_path = ld_path, alpha = alpha, nfolds = nfolds, lambda.min = lambda.min, nlambda = nlambda)

  #run EM to derive epsilon_effects
  epsilon_effect = genee_EM(betas = beta_regularized)

  #run test
  tempresults = genee_loop(betas = beta_regularized, ld = ld, epsilon_effect = epsilon_effect, prior_weight = prior_weight, gene_list = gene_list)

  #return results
  return(tempresults)
}
