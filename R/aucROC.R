### acu.ROC
### This function provide the receiver operating 
### characteristic curve and Area under the Curve of ROC.
### file The connection to the data in the Univariate file.
### pROC roc
### pROC ci.auc
### A matrix to outline the AUC and the confidence region.
aucROC <- function(file) {
  pwdfile = paste(getwd(), "/Univariate/DataTable.csv", sep = "")
  file = pwdfile
  x <- read.csv(file, sep = ",", header = TRUE)
  x.x = x[, 3:ncol(x)]
  rownames(x.x) = x[, 2]
  k = matrix(x[, 1], ncol = 1)
  slink = paste(getwd(), "/PreTable","/slink.csv", sep="")
  slink = read.csv(slink, header=TRUE)
  x.n = cbind(k, x.x)
  sorted = x.n[order(x.n[, 1]), ]
  g = c()
  for (i in 1:nrow(sorted)) {
    if (any(g == sorted[i, 1])) {
      g = g
    }
    else {
      g = matrix(c(g, sorted[i, 1]), ncol = 1)
    }
  }
  NoF = nrow(g)
  dirout.w = paste(getwd(), "/Univariate/ROC", sep = "")
  dir.create(dirout.w)
  
  for (i in 1:NoF) {
    for (j in 1:NoF) {
      
      if (i < j) {
        #pb <- txtProgressBar(min = 1, max = NoF, style = 3)
        
        #Sys.sleep(0.000001)
        
        ni = paste("r.", ExcName(i,slink), ".csv", sep = "")
        nj = paste("r.", ExcName(j,slink), ".csv", sep = "")
        pwdi = paste(getwd(), "/Univariate/Groups/", ni, sep = "")
        pwdj = paste(getwd(), "/Univariate/Groups/", nj, sep = "")
        I = read.csv(pwdi, header = TRUE)
        J = read.csv(pwdj, header = TRUE)
        I = I[, -1]
        J = J[, -1]
        Ilf = matrix(rep(i),nrow(I))
        Jlf = matrix(rep(j),nrow(J))    
        colnames(Ilf) = c("lf")
        colnames(Jlf) = c("lf")
        I = cbind(Ilf,I)
        J = cbind(Jlf,J)
        IJ = rbind(I,J)
        IJM <- as.matrix(IJ[,2:ncol(IJ)])
        outf <-as.factor(IJ[,1])
        message(paste("\n*Group.", ExcName(i,slink), sep = ""),
                " Vs.", paste(" Group.", ExcName(j,slink), sep = ""))
        myROC=function(x,y){
          roc.obj <- pROC::roc(y,x,percent = FALSE)
          auc.ci <- pROC::ci.auc(roc.obj, method = "delong", 
                                 progress = "none", parallel=FALSE)
          rocdata <- c(roc.obj$auc,auc.ci[1],auc.ci[3])
          names(rocdata)=c("AUC","lowAUC","upAUC")
          return(rocdata)
        }
        apply_pb <- function(X, MARGIN, FUN, ...)
        {
          env <- environment()
          pb_Total <- sum(dim(X)[MARGIN])
          counter <- 0
          pb <- txtProgressBar(min = 0, max = pb_Total,
                               style = 3)
          
          wrapper <- function(...)
          {
            curVal <- get("counter", envir = env)
            assign("counter", curVal +1 ,envir= env)
            setTxtProgressBar(get("pb", envir= env),
                              curVal +1)
            FUN(...)
          }
          res <- apply(X, MARGIN, wrapper, ...)
          close(pb)
          res
        }
        
        myroc <- apply_pb(IJM,2,function(x){myROC(x,outf)})
        
        #myroc <- apply(IJM,2,function(x){myROC(x,outf)})
        myroc <- as.data.frame(t(myroc))
        myroc.ij = paste("auc_roc_", ExcName(i,slink), "vs", 
                         ExcName(j,slink), ".csv", sep = "")
        assign(myroc.ij, myroc)
        write.csv(myroc, paste(dirout.w, myroc.ij, sep = "/"))
        #setTxtProgressBar(pb, i)
       
        
      }
    }
  }
}
        