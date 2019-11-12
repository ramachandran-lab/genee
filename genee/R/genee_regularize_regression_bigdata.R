#' Regularized Regression.
#'
#' @param betas A vector of GWAS effect size estimators (OLS betas).
#' @param ld_path Path to the ld matrix memory mapping file (.desc, .bin) including the prefix of the file name. Make sure that .desc and .bin are in the same directory.
#' @param alpha A tuning parameter determining the type of regularization. More specifically, alpha = 0: Ridge Regression; 0 < alpha < 1: Elastic Net; alpha = 1: LASSO.
#' @param nfolds A number indicates how many folds to do for cross validation for regularize regression. The default is 10. If it's fed with 0, then it will perform no cross validation and choose lambda with smallest penalty.
#' @param lambda.min The smallest lambda to use for regularize regression, as a fraction of the largest lambda. Default is 0.0001.
#' @param nlambda The number of lambda values to use. Default is 100.
#' @return A vector of regularized effect size estimators.
#' @export
#' @examples
#' x1 = c(0, 1, 1)
#' x2 = c(0, 1, 2)
#' x3 = c(1, 1, 2)
#' x = cbind(x1 ,x2, x3)
#' ld = cor(x)
#' write.table(ld, file = "myld", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
#' ld_name = "myld"
#' prepare_LD(ld_name)
#' ld_path = "myld"
#' betas = c(0.1, 0.2, 0.1)
#' genee_regularize_regression_bigdata(betas, ld_path, 0.5, 5, 0.001, 50)
genee_regularize_regression_bigdata<-function(betas, ld_path, alpha, nfolds, lambda.min, nlambda){
  # fit a regularized model
  # alpha is the tuning parameter determining which type of the shrinkage is used
  # alpha = 0: Ridge Regression
  # 0 < alpha < 1: Elastic Net
  # alpha = 1: LASSO

  # attach ld matrix
  ld = attach.big.matrix(paste(ld_path, ".desc", sep = ''))
  # Technically, glmnet assumes alpha = 1 for LASSO. However, the LASSO solution can unstable when setting alpha = 1. Thus, we recommend using 0.99 or 0.98 instead.
  if(nfolds == 0){
    fit = biglasso(X = ld, y = betas, alpha = alpha, lambda.min = lambda.min, nlambda = nlambda)

    # use effect size with smallest penalty
    beta_regularized = as.numeric(fit$beta[,ncol(fit$beta)])
  }else{
    cvfit = cv.biglasso(X = ld, y = betas, alpha = alpha,  nfolds = nfolds, lambda.min = lambda.min, nlambda = nlambda)
    # use effect size with smallest error in cross valiation
    beta_regularized = as.numeric(cvfit$fit$beta[,cvfit$min])
  }

  # return results
  return(beta_regularized)
}
