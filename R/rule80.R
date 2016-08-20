#' this function provide the 80 percent rule for the missing value correction
#' @param degree indicate the group factor in the data
#' @param m is the data input
#' @details A variable will be kept if it has a non-zero value for at least 80 percent of each group.
#' @references Smilde A. K.(2005). Fusion of mass spectrometry-based metabolomics data. Anal. Chem. 77, 6729-36.
#' @export
#' @usage  rule80(m,degree)
rule80 <- function(m,degree) {
  dx <- c() 
  for(i in 1:ncol(m)){
    freq <- as.vector(tapply(m[,i] , degree, function(x){sum(x == "0")/length(x)}))
    if(sum(freq > 0.8) > 0) dx <- c(dx , i)
  } 
  if(length(dx) >0) m <- m[,-dx]
  return(m)
}