### Description statistic
### bStatX provide the basic data description for each features. 
### Data descriptions includes mean value, median value, sum, 
### quartile, standard derivatives, etc.
### file The file with the expression information. 
### A matrix for data description
### bStatX(file)
bStatX <- function(file){
    xr <- read.csv(file, sep=",", header=TRUE)
    #xr = file
    xs = xr[,3:ncol(xr)]
    x = cbind(xr[,2],xr[,1],xs)
    x.nn = x
    sorted = x.nn[order(x.nn[, 1]), ]
    g = c()
    for (i in 1:nrow(sorted)) {
          if (any(g == sorted[i, 1])) {
          g = g
          } else {
           g = matrix(c(g, sorted[i, 1]), ncol = 1)
          }
}
dirout.g = paste(getwd(), "/statTarget/statAnalysis/tmp", sep = "")
dir.create(dirout.g)
slink = paste(getwd(), "/statTarget/statAnalysis/PreTable","/slink.csv", sep="")
slink = read.csv(slink, header=TRUE)
for (i in 1:nrow(g)) {
  vuota <- c()
  fin = matrix(rep(NA, ncol(sorted)), nrow = 1)
  for (j in 1:nrow(sorted)) {
    if (sorted[j, 1] == i) {
      vuota <- matrix(sorted[j, ], nrow = 1)
      rownames(vuota) = rownames(sorted)[j]
      fin = rbind(fin, vuota)
    }
  }
  
  nam = paste("r", ExcName(i,slink), sep = ".")
  n = matrix(fin[-1, ], ncol = ncol(sorted))
  n.x = matrix(n[, -1], ncol = ncol(sorted) - 1)
  colnames(n.x) = colnames(x.nn[,2:ncol(x.nn)])
  name = as.matrix(assign(nam, n.x))
  outputfileg = paste("r.", ExcName(i,slink), ".csv", sep = "")
  write.csv(name, paste(dirout.g, outputfileg, sep = "/"), row.names = FALSE)
}
dirout.w = paste(getwd(), "/statTarget/statAnalysis/bstat.Test", sep="")
dir.create(dirout.w)
NoF = nrow(g)
for (i in 1:NoF) {
      ni=paste("r.",ExcName(i,slink),".csv",sep="")
      pwdi = paste(getwd(), "/statTarget/statAnalysis/tmp/", ni, sep="")
      I=read.csv(pwdi, header=TRUE)
      I = I[,-1]
      bS = bStat(I)
      bStat.i=paste("bStat_",ExcName(i,slink), ".csv", sep="")
      assign(bStat.i,bS)
      write.csv(t(bS), paste(dirout.w, bStat.i, sep="/"))
      }
   tmpfile = paste(getwd(), "/tmp/",sep="")
   unlink(tmpfile, recursive=TRUE)
  }










