#' Loop through all gene sets and derive test statistics and p-values for each set.
#'
#' @param betas A vector of effect size estimators (OLS betas or regularized betas).
#' @param ld A matrix of pairwise correlations statistics (i.e., r).
#' @param epsilon_effect The threshold for epsilon-genic effects (i.e., the second largest variance component).
#' @param prior_weight A vector specifying the prior weight for each SNP.
#' @param gene_list A list where each element represents a vector containing the indices of all SNPs in one set (e.g., a set could refer to a gene).
#' @return A list where the first element is a vector of test statistics, and the second element is a vector of p values.
#' @export
#' @examples
#' x1 = c(0, 1, 1)
#' x2 = c(0, 1, 2)
#' x3 = c(1, 1, 2)
#' x = cbind(x1 ,x2, x3)
#' ld = cor(x)
#' betas = c(0.1, 0.2, 0.1)
#' prior_weight = c(1, 1, 1)
#' gene_list=list()
#' gene_list[[1]]=c(1, 2)
#' gene_list[[2]]=c(2, 3)
#' genee_loop(betas, ld, 1e-4, prior_weight, gene_list)
genee_loop <- function(betas, ld, epsilon_effect, prior_weight, gene_list){

  #number of genes
  ngenes = length(gene_list)

  #apply hypothesis testing for each of the gene
  result = unlist(sapply(gene_list, genee_test, ld, betas, epsilon_effect, prior_weight))

  #derive test statistics
  gene_effect_sizes = result[c(1:length(gene_list))*2-1]

  #derive p-values
  gene_p_values = result[c(1:length(gene_list))*2]

  #return results
  return(list(gene_effect_sizes, gene_p_values))
}
