#' @importFrom pROC roc
#' @importFrom pROC ci.auc
auc.ROC <- function(file) {
  pwdfile = paste(getwd(), "/Univariate/DataTable.csv", sep = "")
  file = pwdfile
  x <- read.csv(file, sep = ",", header = TRUE)
  x.x = x[, 3:ncol(x)]
  rownames(x.x) = x[, 2]
  k = matrix(x[, 1], ncol = 1)
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
        ni = paste("r.", i, ".csv", sep = "")
        nj = paste("r.", j, ".csv", sep = "")
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
        myROC=function(x,y){
          roc.obj <- pROC::roc(y,x,percent = FALSE)
          auc.ci <- pROC::ci.auc(roc.obj, method = "bootstrap",boot.n = 500, 
                                 progress = "none")
          rocdata <- c(roc.obj$auc,auc.ci[1],auc.ci[3])
          names(rocdata)=c("AUC","lowAUC","upAUC")
          return(rocdata)
        }
        myroc <- apply(IJM,2,function(x){myROC(x,outf)})
        myroc <- as.data.frame(t(myroc))
        myroc.ij = paste("auc_roc_", i, "vs", j, ".csv", sep = "")
        assign(myroc.ij, myroc)
        write.csv(myroc, paste(dirout.w, myroc.ij, sep = "/"))
      }
    }
  }
}
        