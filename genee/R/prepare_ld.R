#' Preparing LD memory mapping file.
#'
#' @param ld_name The name of LD file.
#' @return Write out the ld memory mapped file (.desc, .bin).
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
prepare_LD<-function(ld_name){
  #prepare ld file using setupX from biglasso
  question <- readline("Warnings: This function will write out a file that has similar size of you data.So make sure you have enough space on the disc. Type yes to continue.")
  if(regexpr(question,'yes')==1){
    ld <- setupX(ld_name, sep = ' ')
  }
}
