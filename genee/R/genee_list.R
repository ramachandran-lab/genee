#' Derive gene list using hg19 genes from UCSC browser.
#'
#' @param glist.hg19 A table containing all the information for hg19 genes from UCSC genome browser.
#' @param all_chr A vector specifying the chromosome for each SNP.
#' @param all_pos A vector specifying the physical location for each SNP.
#' @param upper A number specifying the number of base pairs downstream of the gene to be included.
#' @param lower A number specifying the number of base pairs upstream of the gene to be included.
#' @return A list where the first element is a table containing information of all the genes included in the analysis, and the second element gives a list of indices for each gene.
#' @export
genee_list <- function(glist.hg19, all_chr, all_pos, upper, lower){

  #all genes in the table
  all_genes<-as.character(glist.hg19[,4])

  #find the number of SNPs in each of the gene
  genes_num_snps<-rep(0,length = length(all_genes))

  #making gene_list
  gene_list<-list()
  j=0
  #writing results which contains 7 columns:
  #chr, gene, start, end, nsnps
  geneinfo=matrix(ncol = 5)

  #loop through all the genes in the table
  for (i in 1:length(all_genes)) {
    #current gene
    temp_gene<-all_genes[i]

    #check chromosome (currently, genee doesn't run on sex chromosome)
    if(as.character(glist.hg19[i,1])!="X" && as.character(glist.hg19[i,1])!="Y")
    {
      #get chromosome
      temp_chr=as.numeric(as.character(glist.hg19[i,1]))
      #change gene range according to customer boundary
      start=glist.hg19[i,2]-lower
      end=glist.hg19[i,3]+upper
      #find SNPs included in the ranges
      genes_num_snps[i]=length(which(all_chr==temp_chr & start<=all_pos & end>=all_pos))
      #add gene to the gene_list if there are SNPs in the gene
      if(genes_num_snps[i]!=0 && genes_num_snps[i]!=1){
        j=j+1
        #save the index of SNPs
        gene_list[[j]]=which(all_chr==temp_chr & start<=all_pos & end>=all_pos)
        #temp result for this gene
        tempresults=vector(length = 5)
        #chr
        tempresults[1]=temp_chr
        #gene name
        tempresults[2]=temp_gene
        #start (raw)
        tempresults[3]=start+upper
        #end (raw)
        tempresults[4]=end-lower
        #nsnps
        tempresults[5]=length(gene_list[[j]])
        geneinfo=rbind(geneinfo, tempresults)
      }
    }
  }
  #if no SNPs are foound in the genes, stop
  if(length(gene_list) == 0){
    stop("No gene was found to include any of the target SNPs!")
  }

  #get rid of the first line
  geneinfo = geneinfo[-1,]
  #return results
  return(list(geneinfo, gene_list))
}
