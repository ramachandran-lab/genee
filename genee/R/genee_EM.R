#' EM algorithm to find mixture components and epsilon-genic effects.
#'
#' @param betas A vector of effect size estimators (OLS betas or regularized betas).
#' @return The threshold for epsilon-genic effects (i.e., the second largest variance of all components).
#' @export
#' @examples
#' betas = c(rnorm(1000, 1, 1), rnorm(2000, 2,0))
#' genee_EM(betas)
genee_EM<-function(betas){
  ### EM for normal mixtures ###
  ### This is the function for finding the multiple components of effect sizes and derive epsilon effect sizes for hypothesis testing
  ### genee uses the second largest variance component to construct null hypothesis to test core genes
  # fit a normal mixture model
  # only use non-zero effect sizes to prevent errors during fitting EM
  betas_nonzero = betas[which(betas!=0)]
  
  #if only one element is provided
  if(length(betas_nonzero) <= 1){
    stop("You have no SNP or only 1 SNP has non-zero effect size estimators. Use other regularized method or check data!")
  }

  # fitting univariate normal mixture model
  EM_fit = inference_result_lasso_mclust<-Mclust(betas_nonzero, G = 1:9, modelNames = "V")

  # if only one component returned, then output the variance of the only component as epsilon effect.
  # Otherwise, output the second largest variance.
  if(length(EM_fit$parameters$variance$sigmasq)>1){

    epsilon_effect = sort(EM_fit$parameters$variance$sigmasq, decreasing = TRUE)[2]

  }else{

    epsilon_effect = sort(EM_fit$parameters$variance$sigmasq, decreasing = TRUE)[1]

  }

  #return epsilon effect
  return(epsilon_effect)
}
