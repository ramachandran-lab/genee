# gene-ε

An Association Test for Aggregated Sets of SNP-Level Summary Statistics

## Software Requirements

R (>= 3.3.3), CompQuadForm (>= 1.4.3), glmnet (>= 2.0-16), mclust (>= 5.4.1)

## Installation Instructions

After the ZIP file has been downloaded, type the following command in the terminal:

    R CMD INSTALL genee_0.0.0.9000.tar.gz

Alternatively, one can also install the package from the shell as a [source](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

## Required Function Inputs

* `mydata`: A matrix of summary statistics with four columns. The first column specifies the chromosome for each SNP. The second column specifies the name of each SNP. The third column specifies the physical genomic position for each SNP. The fourth column specifies the GWAS effect size estimator for each SNP.

* `ld`: A matrix of pairwise (Pearson) correlations between all SNPs. Please make sure the SNPs in the LD matrix match/coordinate with the summary statistics table.

## Optional Function Inputs

* `alpha`: A tuning parameter determining the type of regularized regression. Default is alpha = 0.5, which uses the Elastic Net solution. More specifically, alpha = 0: Ridge Regression; 0 < alpha < 1: Elastic Net; alpha = 1: LASSO. Set alpha = -1 to use the original summary statistics with no regularization or shrinkage.

* `upper`: The number of base pairs downstream of the gene to be included. Default is 50000.

* `lower`: The number of base pairs upstream of the gene to be included. Default is 50000.

* `prior_weight`: A vector specifying the prior weight for each SNP. Default assumes equal weights by setting this to be a vector of 1's.

* `gene_list`: A list where each element represents a vector containing the indices of all SNPs in one set (e.g., a set could refer to a gene). If a list is not provided by user, this will be derived using the hg19 gene list from the UCSC browser. In this case, you will have to run the method in the genee package directory in order to load hg19 gene list.

## Tutorial using Simulated Data:

We provide a small simulated example to illustrate how to use gene-ε. Here, we created a simulated dataset using `Simulation.R` file. Under the `Simulated_Data_Example` directory, there are four files:

* `Simulated_LD.RData`

* `Simulated_Summary_Statistics.RData`

* `Simulated_Example.RData`

* `gene_list.RData`

The first two files contain the summary statistics and LD matrix from the simulation. The third file stores all the information about the simulated data. The fourth file is a gene list file. In this case, the function requires that we use this file. To run gene-ε on this simulated dataset, use the following commands:

    ### Load in the LD Matrix ###
    
    load("pathtogenee-master/genee/Simulated_Data_Example/Simulated_LD.RData")
    
    ### Load in the Summary Statistics ###
    load("pathtogenee-master/genee/Simulated_Data_Example/Simulated_Summary_Statistics.RData")

    ### Load in the Gene List ###

    load("pathtogenee-master/genee/Simulated_Data_Example/gene_list.RData")

    ### Run the Analysis ###

    myresult = genee(mydata, ld, gene_list = gene_list)

    ### Load in Info About the Simulation ###
    
    load("pathtogenee-master/genee/Simulated_Data_Example/Simulated_Example.RData")

    ### Create a Manhattan Plot for the Genes ###

    plot(-log10(myresult[,2]), xlab = "gene", ylab = "-log10p", cex = 0.8, pch = 16)

    ### Check for the Causal Genes ###

    points(simu_result$assign_result$all_genes$causal_genes, -log10(myresult[,2])[simu_result$assign_result$all_genes$causal_genes], col = "red",  cex = 0.8, pch = 16)

Note that this analysis uses only default choices for gene-ε. For more details, please see the `Simulation.R` file.

## Tutorial using Real Data:

We provide a small real data example using gene-ε. For this tutorial, we first derive summary statistics for a small set of 2000 SNPs on chromosome 22 using UK Biobank (self-identified as British) individual-level genotype data and the `Body Height` trait. Next, we used the CEU (Western European Ancestry) and GBR (British in England and Scotland) population genotype data from the 1000 Genomes Project to derive a corresponding LD matrix. To run gene-ε on this dataset, use the following commands:

    ### Set the Working Directory ###
    
    setwd("pathtogenee-master/genee/") 
    
    #NOTE: Make sure to run this within the genee directory to load the hg19 gene list file!

    ### Load in the LD Matrix ###

    load("pathtogenee-master/genee/Real_Data_Example/Real_LD.RData")

    ### Load in the Summary Statistics ###
    load("pathtogenee-master/genee/Real_Data_Example/Real_Summary_Statistics.RData")
    
    ### Run the Analysis ###

    myresult = genee(mydata, ld)

Again, this analysis uses the default choices for gene-ε.

## Additional Notes and Thoughts:

* As you may notice, gene-ε does not require a gene list file. It can derive region based p-values for  SNPs contained within boundaries that users set up themselves. Therefore, users can provide any list that allows them to test arbitrary sets of SNPs.

* gene-ε requires an LD matrix as an input. Hence, if users have large genome-wide dataset (i.e. a million SNPs), it may be a struggle for R to handle such a large LD matrix. Thus, in these cases, we recommend considering to run the method on a local or *cis* basis (e.g., chromosome by chromosome).

## Questions and Feedback: 
For questions or concerns with gene-ε, please contact [Wei Cheng](mailto:wei_cheng1@brown.edu).

We appreciate any feedback you may have with our software and/or instructions.
