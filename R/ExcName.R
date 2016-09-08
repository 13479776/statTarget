ExcName <- function(i,slink){
  e <- grep(i,slink[,2])
  slink[e[1],1]
}