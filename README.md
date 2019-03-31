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

1.mydata : A matrix of four columns. The first column specifies the chromosome for each SNP. The second column specifies the name for each SNP. The third column specifies the physical position for each SNP. The fourth column specifies the GWAS effect size estimator for each SNP.

2.ld: A matrix of pairwise correlations statistics (i.e., r). Note: Please make sure the SNPs in LD matrix and in the summary statistics table match.

# Optional Input:

3.alpha: A tuning parameter determining the type of regularized regression. Default is alpha = 0.5 which uses the Elastic Net solutitableon. More specifically, alpha = 0: Ridge Regression; 0 < alpha < 1: Elastic Net; alpha = 1: LASSO. Set alpha = -1 to use the original summary statistics with no regularization or shrinkage.

4.upper: The number of base pairs downstream of the gene to be included. Default is 50000.

5.lower: The number of base pairs upstream of the gene to be included. Default is 50000.

6.prior_weight:	A vector specifying the prior weight for each SNP. Default assumes equal weights by defining prior_weight to be a vector containing 1's.

7.gene_list:	A list where each element represents a vector containing the indices of all SNPs in one set (e.g., a set could refer to a gene). If it's not provided by user, this will be derived using the hg19 gene list from UCSC browser. If this parameter is not provided, you will have to run the method in the genee package directory in order to load hg19 gene list.

# Example:
We provide two small example datasets to illustrate how to use genee.

1.Real example data: Under the example_data_real directory, there are two files. "example_ld_real.RData" and "example_summary_statistics_real.RData". This dataset is created using real data. We first derive summary statistics for a small set of SNPs (2000 SNPs) on chromosome 22 for Height using UKBiobank (self-identified as British) genotype and phenotype data. Then we use 1000 genome CEU (Western European Ancestry) and GBR (British in England and Scotland) populatio genotype data to derive corresponding LD matrix. To run genee on this dataset, you will need the following commands.

setwd("pathtogenee-master/genee/") # Make sure to run this within the genee directory to load the hg19 gene list file!"

load("pathtogenee-master/genee/example_data_real/example_ld_real.RData")#This will give you ld matrix

load("pathtogenee-master/genee/example_data_real/example_summary_statistics_real.RData")#This will give you mydata

myresult = genee(mydata, ld)#Using default parameters, alpha = 0.5, upper = 50000, lower = 50000, prior_weight = vector of 1's.



2.Simulate example data: We also create simulate data set for you. The data is generated using Simulation.R (We also encourage users to generate simulate data by their owns using this simulation strategy or their own simulations). Under the example_data_simulate directory, there are four files. "example_ld_simulate.RData", "example_summary_statistics_simulate.RData", "example_simulate.RData" and "gene_list.RData". The first two files are basically similar to real example data. The third file is a list stores all the information about this simulation data (For more details, please see Simulation.R file). The fourth file is a gene_list file that will also be required in this case. To run genee on this dataset, you will need the following commands.

load("pathtogenee-master/genee/example_data_simulate/example_ld_simulate.RData")

load("pathtogenee-master/genee/example_data_simulate/example_summary_statistics_simulate.RData")

load("pathtogenee-master/genee/example_data_simulate/gene_list.RData")

myresult = genee(mydata, ld, gene_list = gene_list)#Using default parameters, alpha = 0.5, upper = 50000, lower = 50000, prior_weight = vector of 1's.

load("pathtogenee-master/genee/example_data_simulate/example_simulate.RData")#load info about the simulation

plot(-log10(myresult[,2]), xlab = "gene", ylab = "-log10p", cex = 0.8, pch = 16) # Manhattan plot

points(simu_result$assign_result$all_genes$causal_genes, -log10(myresult[,2])[simu_result$assign_result$all_genes$causal_genes], col = "red",  cex = 0.8, pch = 16) # check real causal genes



# Note:
1. As you may notice, genee doesn't require a gene_list file. It can derive genes containing the SNPs according to the boundary that users set up using hg19 gene list. Therefore, to use hg19 gene list, users will have to run genee wihthin genee directory to load the glist-hg19 file. Meanwhile, users can also provide the gene_list file by themselves for testing arbitrary sets of SNPs.

2. genee will require a LD matrix as input. Hence if users have large genome wide dataset (i.e. a million SNPs), it will take a long time for R to handle such a huge matrix. Thus, we recommend to run the method by chromosome in this case.
