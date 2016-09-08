RsdCal = function(file,batch = TRUE, DistPattern = TRUE, output= TRUE){
  #xr <- read.csv(file, sep=",", header=TRUE)
  xr = file
  #xr = xr[,-1]
  xs = xr[,3:ncol(xr)]
  if(batch){
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
    dirout.g = paste(getwd(), "/tmp", sep = "")
    dir.create(dirout.g)
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
      nam = paste("r", i, sep = ".")
      n = matrix(fin[-1, ], ncol = ncol(sorted))
      n.x = matrix(n[, -1], ncol = ncol(sorted) - 1)
      colnames(n.x) = colnames(x.nn[,2:ncol(x.nn)])
      name = as.matrix(assign(nam, n.x))
      outputfileg = paste("r.", i, ".csv", sep = "")
      write.csv(name, paste(dirout.g, outputfileg, sep = "/"), 
                row.names = FALSE)
    }
    dirout.w = paste(getwd(), "/statTarget/shiftCor/RSDresult", sep="")
    dir.create(dirout.w)
    NoF = nrow(g)
    nrow = dim(xr)[2] -2
    batchQC <- matrix(rep(NA),ncol = NoF,nrow = nrow)
    for (i in 1:NoF) {
      ni=paste("r.",i,".csv",sep="")
      pwdi = paste(getwd(), "/tmp/", ni, sep="")
      I=read.csv(pwdi, header=TRUE)
      I = I[,-1]
      bS = t(bStatCor(I)[nrow(bStatCor(I)),])
      bS = as.matrix(bS)
      batchQC[,i] = bS
    }
    bSall = t(bStatCor(xs)[nrow(bStatCor(xs)),])
    batchQC <- cbind(batchQC,bSall)
  rownames(batchQC) <- NULL
  colnames(batchQC) <- NULL
  ###########RSDdist##########################
  #message("Summary QC distribution of the CV 
  #in each batch and all batch (last one):")
  distout <- matrix(rep(NA),nrow = NoF+1,ncol = 20)
  for(i in 1:(NoF+1)) {
        valueB <- batchQC[,i]
        Batchdist <- sapply(seq(0,95,5),function(x){
          sum(valueB*100 <=x+5,na.rm = TRUE)*100/length(valueB)})
        distout[i,] <- matrix(Batchdist,1)
  }
  colnames(distout) <- paste("CV","<",seq(5,100,5),"%",sep = "")
  rownames(distout) <- paste("Batch","_",seq(1,nrow(distout),1),sep = "")
  rownames(distout)[nrow(distout)] <- "Total"
  
  if(DistPattern){
    print(distout)
  }
  
  
  #return(distout)
  
  if(output){
    RSDdist_CV=paste("RSDdist_QC_stat",".csv", sep="")
    write.csv(distout, paste(dirout.w, RSDdist_CV, sep="/"))
    RSD_all=paste("RSD_all",".csv", sep="")
    write.csv(bSall, paste(dirout.w, RSD_all, sep="/"))
  }
  tmpfile = paste(getwd(), "/tmp/",sep="")
  unlink(tmpfile, recursive=TRUE)
  } else {bSall = t(bStatCor(xs)[nrow(bStatCor(xs)),])}
  
  return(bSall)
  
}
