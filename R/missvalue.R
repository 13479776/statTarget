#' Miss value will be replaced with median value.
#' @param degree is the group factor
#' @param mv is th data file
#' @usage  missvalue(mv, degree)
missvalue <- function(mv,degree) {
  for (i in 1:dim(mv)[1]){  
    for(j in 2:dim(mv)[2]){
      if(mv[i,j] == "0" ){
        
        mv[i,j] <- tapply(mv[,j],degree,median,na=TRUE)[degree[i]]
        
      }
    }
  } 
  return(mv)  
}
