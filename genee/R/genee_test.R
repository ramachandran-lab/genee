#' Hypothesis Testing for One Set of SNPs
#'
#' @param gene A vector containing the indices of all the SNPs in the set (e.g., a set may refer to a gene).
#' @param ld A matrix of pairwise correlations statistics (i.e., r).
#' @param betas A vector of effect size estimators (OLS betas or regularized betas).
#' @param epsilon_effect The threshold for epsilon-genic effects (i.e., the second largest varaince of all components).
#' @param prior_weight A vector specifying the prior weight for each SNP.
#' @return A vector containing test statistic, variance and p-value for the set.
#' @export
#' @examples
#' gene = c(1,2)
#' x1 = c(0, 1, 1)
#' x2 = c(0, 1, 2)
#' x3 = c(1, 1, 2)
#' x = cbind(x1 ,x2, x3)
#' ld = cor(x)
#' betas = c(0.1, 0.2, 0.1)
#' prior_weight = c(1, 1, 1)
#' genee_test(gene, ld, betas, 1e-4, prior_weight)
genee_test<-function(gene, ld, betas, epsilon_effect, prior_weight){

  #local ld
  ld_g = as.matrix(ld[gene, gene])

  #number of SNPs in the gene
  nsnp_g = length(gene)

  #define weight matrix
  weight_matrix <- diag(x = prior_weight[gene], nsnp_g, nsnp_g)

  #derive eigenvalues
  temp_e <- eigen((ld_g*epsilon_effect)%*%weight_matrix)
  e_values = temp_e$values
  

  #compute test statistics
  test_statsics = betas[gene]%*%weight_matrix%*%betas[gene]

  #compute test statistics variance
  t_var = sum(diag((ld_g*epsilon_effect)%*%(ld_g*epsilon_effect)))

  #Using imhof to derive p-values
  options(warn=1)
  ans=imhof(test_statsics, lambda = e_values,h = rep(1, length(e_values)), delta = rep(0,length(e_values)), epsabs = 10^(-16), epsrel = 10^(-16), limit = 1000)

  #prevent negative p-values
  p_value_g <-ans$Qq
  if(p_value_g <= 0){
    p_value_g = 10^-20
  }

  #return test statistics and p-values
  return(c(test_statsics, t_var, p_value_g))
}
