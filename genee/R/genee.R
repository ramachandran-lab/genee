#' genee: An Association Test for Aggregated Sets of SNP-Level Summary Statistics.
#'
#' @param mydata A matrix of four columns. The first column specifies the chromosome for each SNP. The second column specifies the name for each SNP. The third column specifies the physical position for each SNP. The fourth column specifies the GWAS effect size estimator for each SNP.
#' @param ld A matrix of pairwise correlations statistics (i.e., r).
#' @param alpha A tuning parameter determining the type of regularized regression. Default is alpha = 0.5 which uses the Elastic Net solution. More specifically, alpha = 0: Ridge Regression; 0 < alpha < 1: Elastic Net; alpha = 1: LASSO. Set alpha = -1 to use the original summary statistics with no regularization or shrinkage.
#' @param upper The number of base pairs downstream of the gene to be included. Default is 50000.
#' @param lower The number of base pairs upstream of the gene to be included. Default is 50000.
#' @param prior_weight A vector specifying the prior weight for each SNP. Default assumes equal weights by defining prior_weight to be a vector containing 1's.
#' @param gene_list A list where each element represents a vector containing the indices of all SNPs in one set (e.g., a set could refer to a gene). If it's not provided by user, this will be derived using the hg19 gene list from UCSC browser.
#' @param nfolds A number indicates how many folds to do for cross validation for regularize regression. The default is 10. If it's fed with 0, then it will perform no cross validation and choose lambda with smallest penalty.
#' @return If a gene list is not provided, genee will return a matrix containing the following 7 columns: 1. chromosome number (chr); 2. gene name; 3. start position; 4. end position; 5. number of SNPs in the gene/set (nsnp); 6. gene/set test statistics (test_q); 7.gene/set test statistics variance and 8. gene/set p-value (pval). If a gene list is provided, genee will return a matrix with three columns: 1. gene/set test statistics (test_q);  2.gene/set test statistics variance and 3. gene/set p-value (pval).
#' @export
#' @examples
#' x1 = c(0, 1, 1)
#' x2 = c(0, 1, 2)
#' x3 = c(1, 1, 2)
#' x = cbind(x1 ,x2, x3)
#' ld = cor(x)
#' prior_weight = c(1, 1, 1)
#' gene_list=list()
#' gene_list[[1]]=c(1, 2)
#' gene_list[[2]]=c(2, 3)
#' chr = c(1, 1, 1)
#' snpnames = c('rs1', 'rs2', 'rs3')
#' pos = c(100, 200, 300)
#' betas = c(0.1, 0.2, 0.1)
#' mydata = matrix(nrow = 3, ncol = 4)
#' mydata[,1] = chr
#' mydata[,2] = snpnames
#' mydata[,3] = pos
#' mydata[,4] = betas
#' genee(mydata, ld, alpha =-1, gene_list = gene_list)
genee <- function(mydata, ld, alpha, upper, lower, prior_weight, gene_list, nfolds){

  #get chromsomes
  all_chr=as.numeric(mydata[,1])
  #get SNP names
  all_snps=mydata[,2]
  #get SNP positions
  all_pos=as.numeric(mydata[,3])
  #get GWAS betas
  betas_ols=as.numeric(mydata[,4])

  #make sure ld is a matrix
  ld=as.matrix(ld)

  #if missing alpha, set default to be 0.5
  if(missing(alpha)){
    alpha = 0.5
  }

  #if missing upper, set default to be 50000
  if(missing(upper)){
    upper = 50000
  }

  #if missing lower, set default to be 50000
  if(missing(lower)){
    lower = 50000
  }

  #if missing prior weight, set default to be 1 for each SNP
  if(missing(prior_weight)){
    prior_weight = rep(1, length(betas_ols))
  }

  #if gene_list is not provided by user, run genee_list
  #get gene_list
  if(missing(gene_list)){
    #load gene tables
    load("data/glist.hg19.rda")

    #derive gene list and gene info
    temp = genee_list(glist.hg19 = glist.hg19, all_chr = all_chr, all_pos = all_pos, upper = upper, lower = lower)
    gene_info = temp[[1]]
    gene_list = temp[[2]]
  }else{
    #for given gene list, only output test statistics and pvals
    gene_info = NA
  }

  #if missing nfolds, set default to be 10
  if(missing(nfolds)){
    nfolds = 10
  }

  #if alpha == -1, use raw GWAS betas
  #else using regularized version
  if(alpha == -1){
    tempresults = genee_ols(betas_ols = betas_ols, ld = ld, prior_weight = prior_weight, gene_list = gene_list)
  }else{
    tempresults = genee_regularize(betas_ols = betas_ols, alpha = alpha, ld = ld, prior_weight = prior_weight, gene_list = gene_list, nfolds = nfolds)
  }

  #test statistics
  test_statistics = tempresults[[1]]

  #pvals for genes
  tvar = tempresults[[2]]

  #pvals for genes
  pvals = tempresults[[3]]

  #check whether have gene_info
  if(all(is.na(gene_info))){
    #only output test statistics and pvals
    test_q = test_statistics
    q_var = tvar
    pval = pvals
    final_results = data.frame(test_q, q_var, pval)
  }else{
    chr = as.numeric(gene_info[,1])
    gene = as.character(gene_info[,2])
    start = as.numeric(gene_info[,3])
    end = as.numeric(gene_info[,4])
    nsnp = as.numeric(gene_info[,5])
    test_q = test_statistics
    q_var = tvar
    pval = pvals
    final_results = data.frame(chr, gene, start, end, nsnp, test_q, q_var, pval)
  }

  #return results
  return(final_results)
}

