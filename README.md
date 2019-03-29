# genee
An Association Test for Aggregated Sets of SNP-Level Summary Statistics

# Requirements:
R (>= 3.3.3), CompQuadForm (>= 1.4.3), glmnet (>= 2.0-16), mclust (>= 5.4.1)

# Questions: 
Contact wei_cheng1@brown.edu

# Installation Instructions:
After download the ZIP file, type in the following commands in the terminal:

R CMD INSTALL genee_0.0.0.9000.tar.gz

# Required Input:

1.Summary Statistics table: A matrix of four columns. The first column specifies the chromosome for each SNP. The second column specifies the name for each SNP. The third column specifies the physical position for each SNP. The fourth column specifies the GWAS effect size estimator for each SNP.

2.LD matrix: A matrix of pairwise correlations statistics (i.e., r). Note: Please make sure the SNPs in LD matrix and in the summary statistics table match.

# Optional Input:

3.alpha: A tuning parameter determining the type of regularized regression. Default is alpha = 0.5 which uses the Elastic Net solution. More specifically, alpha = 0: Ridge Regression; 0 < alpha < 1: Elastic Net; alpha = 1: LASSO. Set alpha = -1 to use the original summary statistics with no regularization or shrinkage.

4.upper: The number of base pairs downstream of the gene to be included. Default is 50000.

5.lower: The number of base pairs upstream of the gene to be included. Default is 50000.

6.prior_weight:	A vector specifying the prior weight for each SNP. Default assumes equal weights by defining prior_weight to be a vector containing 1's.

7.gene_list:	A list where each element represents a vector containing the indices of all SNPs in one set (e.g., a set could refer to a gene). If it's not provided by user, this will be derived using the hg19 gene list from UCSC browser.

# Example code:

genee(summary data, ld)

genee(summary data, ld, -1)

genee(summary data, ld, upper=10000, lower=1000)

genee(summary data, ld, 0.98, gene_list = my_list)

genee(summary data, ld, alpha = 0, prior_weight = my_prior)

# Note:
1. As you may notice, genee doesn't require a gene_list file. It can derive genes containing the SNPs according to the boundary that users set up using hg19 gene list. However, users can provide the gene_list file for testing arbitrary sets of SNPs.

2. genee will require a LD matrix as input. Hence if users have large genome wide dataset (i.e. a million SNPs), it will take a long time for R to handle such a huge matrix. Thus, we recommend to run the method by chromosome in this case.
