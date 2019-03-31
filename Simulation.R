#Simulation for Genee

### Clear Environment ###
rm(list=ls())

####random assign causal genes and causal snps only in genes with unique snps
assign_random<-function(gene_list, p_gene, p_snp,ncausal_intergenic){
  #get the number of SNPs
  nsnp=dim(X)[2]
  #get the number of genes
  ngene=length(gene_list)-1
  #get the number of causal genes
  ncausal_gene=round(ngene*p_gene+0.5)
  #sample genes
  causal_genes=sample(c(1:ngene), ncausal_gene)#id of causal genes
  #all genes id
  genes.id=c(1:ngene)
  #non-causal genes
  noncausal_genes=genes.id[-causal_genes]
  #total number of causal SNPs in genes
  ncausal_snps_in_genes=round(nsnp*p_snp+0.5)
  #id of all snps 
  snps.id=c(1:nsnp)
  #id of SNPs in causal genes
  snps_in_causal_genes=c()
  for (j in causal_genes) {
    snps_in_causal_genes=c(snps_in_causal_genes, gene_list[[j]])
  }
  #assert that total SNPs in genes are more than total number of causal SNPs
  stopifnot(length(snps_in_causal_genes)>ncausal_snps_in_genes)
  #sample causal SNPs in genes
  causal_snps_in_genes<-sample(snps_in_causal_genes, ncausal_snps_in_genes)
  #snps in the intergenic region
  snps_intergenic=snps.id[c((max(gene_list[[ngene]])+1):nsnp)]#id of intergenic snps
  #assert that the number of intergenic snps are more than causal snps intergenic region
  stopifnot(length(snps_intergenic)>ncausal_intergenic)
  #sample intergenic causal SNPs
  causal_snps_intergenic<-sample(snps_intergenic, ncausal_intergenic)#id of causal snps in intergenic
  #all causal snps
  causal_snps=sort(unique(c(causal_snps_in_genes, causal_snps_intergenic)))
  #all non-causal snps
  noncausal_snps=snps.id[-causal_snps]
  #some gene might not be sampled
  ncausal_gene=rep(0,(length(gene_list)-1))
  for (i in 1:(length(gene_list)-1)) {
    ncausal_gene[i]=length(intersect(gene_list[[i]], causal_snps))
  }
  #causal genes id
  causal_genes=which(ncausal_gene!=0)
  #non-causal genes id
  noncausal_genes=which(ncausal_gene==0)
  #store results
  all_genes=list()
  all_genes$causal_genes=causal_genes
  all_genes$noncausal_genes=noncausal_genes
  all_snps=list()
  all_snps$causal_snps=causal_snps
  all_snps$causal_snps_in_genes=causal_snps_in_genes
  all_snps$causal_snps_intergenic=causal_snps_intergenic
  all_snps$noncausal_snps=noncausal_snps
  assign_result=list()
  assign_result$all_genes=all_genes
  assign_result$all_snps=all_snps
  return(assign_result)
}


########################################################
########################################################
###############phenotype simulation#####################

pheno_simu<-function(eta_total, pve, X, assign_result, PCs){
  ind=dim(X)[1]
  nsnp=dim(X)[2]
  Xmarginal=X[,assign_result$all_snps$causal_snps]
  Xnoise=X[,assign_result$all_snps$noncausal_snps]
  #marginal effects
  beta_marginal=rnorm(dim(Xmarginal)[2], 0, sd=1)
  y_marginal=Xmarginal%*%beta_marginal
  beta_marginal=beta_marginal*as.numeric(sqrt(pve*(1-eta_total)/var(y_marginal)))
  y_marginal=Xmarginal%*%beta_marginal
  #noisy efffects if eta_total!=0
  beta_noise=rnorm(dim(Xnoise)[2], 0, sd=1)
  y_noise=Xnoise%*%beta_noise
  if(eta_total==0){
    y_noise=y_noise*0
    beta_noise=beta_noise*0
  }else{
    beta_noise=beta_noise*as.numeric(sqrt(pve*eta_total/var(y_noise)))
    y_noise=Xnoise%*%beta_noise 
  }
  #error effects
  y_err=rnorm(ind)
  y_err=y_err*as.numeric(sqrt((1-pve)/var(y_err)))
  #y without structure
  y=y_marginal+y_noise+y_err
  true_betas=vector(length=nsnp)
  true_betas[assign_result$all_snps$causal_snps]=beta_marginal
  true_betas[assign_result$all_snps$noncausal_snps]=beta_noise
  pheno_result=list()
  pheno_result$beta_marginal=beta_marginal
  pheno_result$beta_noise=beta_noise
  pheno_result$true_betas=true_betas
  pheno_result$y=y
  pheno_result$y_err=y_err
  pheno_result$y_marginal=y_marginal
  pheno_result$y_noise=y_noise
  ###pcs
  ### Define the effects of the PCs (we set the structure to explain 10 percent of varaince)###
  beta_pc=rnorm(dim(PCs)[2])
  y_pcs=PCs%*%beta_pc
  beta_pc=beta_pc*as.numeric(sqrt(0.1/var(y_pcs)))
  y_pcs=PCs%*%beta_pc
  #error effects
  y_err_wPCs=rnorm(ind)
  y_err_wPCs=y_err_wPCs*as.numeric(sqrt((1-pve-0.1)/var(y_err_wPCs)))
  y_wPCs=y_marginal+y_noise+y_err_wPCs+y_pcs
  pheno_result$wPCs$beta_marginal_wPCs=beta_marginal
  pheno_result$wPCs$beta_noise_wPCs=beta_noise
  pheno_result$wPCs$true_betas_wPCs=true_betas
  pheno_result$wPCs$y_wPCs=y_wPCs
  pheno_result$wPCs$y_pcs=y_pcs
  pheno_result$wPCs$y_err_wPCs=y_err
  pheno_result$wPCs$y_marginal_wPCs=y_marginal
  pheno_result$wPCs$y_noise_wPCs=y_noise
  return(pheno_result)
}



###Linear Regression ###
summary_stats<-function(X, y, PCs){
  num_snps=dim(X)[2]
  pval=vector(length = num_snps)
  beta_est=vector(length = num_snps)
  beta_se=vector(length = num_snps)
  for (i in 1:num_snps) {
    relation=lm(y~X[,i])
    beta_est[i]=summary(relation)$coefficient[2,1]
    beta_se[i]=summary(relation)$coefficient[2,2]
    pval[i]=summary(relation)$coefficient[2,4]
  }
  pval_wPCs=vector(length = num_snps)
  beta_est_wPCs=vector(length = num_snps)
  beta_se_wPCs=vector(length = num_snps)
  for (i in 1:num_snps) {
    relation_wPCs=lm(y~X[,i]+PCs)
    beta_est_wPCs[i]=summary(relation_wPCs)$coefficient[2,1]
    beta_se_wPCs[i]=summary(relation_wPCs)$coefficient[2,2]
    pval_wPCs[i]=summary(relation_wPCs)$coefficient[2,4]
  }
  gwas_result=list()
  gwas_result$beta_est=beta_est
  gwas_result$beta_se=beta_se
  gwas_result$pval=pval
  gwas_result$beta_est_wPCs=beta_est_wPCs
  gwas_result$beta_se_wPCs=beta_se_wPCs
  gwas_result$pval_wPCs=pval_wPCs
  return(gwas_result)
}


########################################################
###############My simulation#######################
my_simulator<-function(X, p_snp, p_gene, ncausal_intergenic, gene_list, eta_total, pve, PCs){
  assign_result=assign_random(gene_list, p_gene, p_snp,ncausal_intergenic)
  pheno_result=pheno_simu(eta_total, pve, X, assign_result, PCs)
  gwas_result=summary_stats(X, pheno_result$y, PCs)
  gwas_result_wPCs_pheno=summary_stats(X, pheno_result$wPCs$y_wPCs, PCs)
  ########################################################
  ###############Summary Statistics#######################
  simu_result=list()
  simu_result$assign_result=assign_result
  simu_result$pheno_result=pheno_result
  simu_result$gwas_result=gwas_result
  simu_result$gwas_result_wPCs_pheno=gwas_result_wPCs_pheno
  return(simu_result)
}

### Simulate gene list
# (1) ngene = number of genes
# (2) min_gsize = minimum number of SNPs in a gene
# (3) max_gsize = maximum number of SNPs in a genes
# (4) nsnp_intergenic = number of SNPs in intergenic region
ngene=100; min_gsize=5; max_gsize=20; nsnp_intergenic=400

### Save the index of SNPs for each gene
gene_list=c()
counter=0
for(i in 1:ngene){
  #simulate the size for a gene
  temp_gsize=sample(c(min_gsize:max_gsize),1)
  #store the index
  gene_list[[i]]=c((counter+1):(counter+temp_gsize))
  counter=counter+temp_gsize
}

### Save the index of SNPs for the inter-genic SNPs
gene_list[[ngene+1]] = c((max(gene_list[[ngene]])+1):(max(gene_list[[ngene]])+nsnp_intergenic))

### Simulate genotypes:
# (1) ind= # number of samples.
# (2) nsnp = number of SNPs or variants.

ind = 3e3; nsnp = max(gene_list[[ngene+1]]);

### Simulate the genotypes such that all variants have minor allele frequency (MAF) > 0.05 ###
# NOTE: As in the paper, we center and scale each genotypic vector such that every SNP has mean 0 and standard deviation 1.
maf <- 0.05 + 0.45*runif(nsnp)
X   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
X   <- matrix(as.double(X),ind,nsnp,byrow = TRUE)
Xoriginal=X
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

### PCs of X
fit_pcs=prcomp(X)
PCs=fit_pcs$x[,c(1:5)]
rm(fit_pcs)

### simulate phenotype
# (1) p_snp= percent of causal SNPs in genes.
# (2) p_gene = percent of causal genes. Note: You have to play around with first two parameters such that there are enough SNPs in causal genes to meet the percent of causal SNPs.
# (3) ncausal_intergenic = number of causal SNPs in intergenic region.
# (4) eta_total = the percent of total variance explained by all non-causal SNPs/noise effect. Default to be 0.
# (5) pve = narrow sense heritability/proportion of varaince explained by additive genetics effects. 
p_snp=0.01; p_gene=0.02; ncausal_intergenic=10; eta_total=0; pve=0.6
#generate data, if error appears, play around with the parameters above and rerun it multiple times can fix it.
simu_result=my_simulator(X, p_snp, p_gene, ncausal_intergenic, gene_list, eta_total, pve, PCs)


### compute LD 
ld=cor(X)


### prepare summary statistics for genee format
example_input=matrix(nrow = nsnp, ncol = 4)
#chr
example_input[,1]=rep(1,nsnp)
#name
for(i in 1:nsnp){
  example_input[i,2]=paste("rs",i, sep="")
}

#pos
example_input[,3]=c(1:nsnp)
#GWAS effect sizes(OLS)
example_input[,4]=simu_result$gwas_result$beta_est_wPCs
mydata = example_input


### Save example
save(simu_result, file = "Simulated_Example.RData")
save(mydata, file = "Simulated_Summary_Statistics.RData")
save(ld, file = "Simulated_LD.RData")
#save gene list (deleting the intergenic region)
gene_list[[ngene+1]] = NULL
save(gene_list, file = "gene_list.RData")
