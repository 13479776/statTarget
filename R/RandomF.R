#' RandomF provide the random forest function.
#' @param nvarRF shows the number of variables in Gini plot of Randomforest model (=< 100). 
#' @param file is the data
#' @description Multi result was outputed, such as MDS plot, vriable Gini plot and list results
#' @usage RandomF(file, nvarRF)
#' @export
RandomF <- function(file,nvarRF) {
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
  dirout.w = paste(getwd(), "/Univariate/RForest", sep = "")
  dirout.RF = paste(getwd(), "/Univariate/RForest/RF_Plots/", sep="")
  dir.create(dirout.w)
  dir.create(dirout.RF)
  
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
        #fin = ncol(sorted) - 1
        #we <- matrix(rep(NA, fin))
        #odds.radio..........................................................................................................................
        IJM <- as.matrix(IJ[,2:ncol(IJ)])
        outf <-as.factor(IJ[,1])
         
        #data(iris)
        #y <- as.factor(ifelse(iris$Species == "setosa" |
        #                        iris$Species == "virginica", 1, 0) )
        #xdata <- iris[,1:4]
        #yf <- c(1,1,2,2)
        #y <- as.factor(yf)
        nrep = 100
        mods <- rlply(nrep, randomForest(IJM, outf, ntree=500))
        dat <- matrix(rep(NA),nrep,length(IJM[1,]))
        #dat1 <- matrix(rep(NA),1,)
        dat1 <-c()
        for (k in 1:nrep) {
          for (h in 1:length(IJM[1,])){
            p <- function(x,y) { dat = as.list(importance(mods[[x]]))[[y]]}
            dat[k,h] <- p(k,h)
            dat1[k] <- median(mods[[k]][[4]][,1])
          }
          colnames(dat) = colnames(IJM)
          dat1 = matrix(dat1)
          colnames(dat1) = c("error.rate")
        }
        #write.csv(dat,"Gini_RF.csv")
        #write.csv(dat1,"error.rate.csv")
        Gini.ij = paste("Gini_", i, "vs", j, ".csv", sep = "")
        assign(Gini.ij, dat)
        write.csv(dat, paste(dirout.w, Gini.ij, sep = "/"))
        error.rate.ij = paste("error.rate_", i, "vs", j, ".csv", sep = "")
        assign(error.rate.ij, dat1)
        write.csv(dat1, paste(dirout.w, error.rate.ij, sep = "/"))
        
        ########### Randomforest MDS,boxplot ######
        
        rfmods <- randomForest(IJM, outf, ntree=500,proximity=TRUE)
        MDS = paste(dirout.RF, "MDSPlot_", i, "vs", j, ".pdf", sep="")
        pdf(MDS)
        rfplot <- MDSplot(rfmods, outf,
        bg=rainbow(length(levels(as.factor(outf))))[unclass=outf], 
        pch=c(21:25)[unclass=outf], 
        palette=rep("black",5), 
        xlab="Dimension 1", 
        ylab="Dimension 2", 
        cex=0.8, main=paste("MDS Plot ", i, " vs ", j, sep=""))
        legend(
          "bottomright",
          levels(outf), 
          pch=c(21:25,21:25), 
          col=rep("black",length(levels(outf))), 
          pt.bg=rainbow(length(levels(outf))), 
          cex=0.6)
        #rfplot <- MDSplot(rfmods, outf,palette=rep(1, 2), pch=as.numeric(outf),col = palette[as.numeric(outf)])
        #plot(rfplot$points,pch=20,bg="black",main=c("Random Forest"),col=c("red","green"))
        #text(Dim1, Dim2, rownames(euro.mds), cex=0.8, col="red")
        dev.off()
        RFscore.ij = paste("RFscore_", i, "vs", j, ".csv", sep = "")
        assign(RFscore.ij, rfplot$points)
        write.csv(rfplot$points, paste(dirout.w, RFscore.ij, sep = "/"))
        Giniplot = paste(dirout.RF, "Giniplot_", i, "vs", j, ".pdf", sep="")
        pdf(Giniplot)
        orderindex <- apply(dat,2,median)
        datindex <- cbind(orderindex,t(dat))
        sort.dat <- datindex[order(datindex[,1], na.last=NA, decreasing = T),]
        sort.dat <- t(sort.dat[,-1]) 
        sort.range <- sort.dat[,1:nvarRF]
        sort.datFilter <- flipdim_sT(sort.range, 2)
        par(mar=c(3, 10, 3, 10) + 0.1,mgp=c(1.5,0.5,0))
        colfunc <- colorRampPalette(c("white", "green"))
        #colfunc <- colorRampPalette(colors = brewer.pal(9,"reds"))
        boxplot(sort.datFilter, main = paste("Gini Plot ", i, " vs ", j, sep=""),
                notch = FALSE, col =colfunc(nvarRF), 
                cex =0.1, las = 1, xlab= "Mean Decrease Gini",ylab="Variables", horizontal = TRUE,
                pars = list(boxwex = 0.6, staplewex = 0.3, outwex = 0.3),lwd=0.3,
                cex.lab=0.8, cex.axis =0.4,medlwd=0.3, xaxt='n',yaxt="n")
        axis(2,labels= c(colnames(sort.datFilter)), at = 1:nvarRF, 
             cex = 0.3, tck= -0.015,cex.axis=0.3,las=1)
        axis(1,cex = 0.5, font= 0.05, tck= -0.015,cex.axis=0.3)
        grid(NULL,NA, lwd = 1,col = "gray")
        #axis(side=1, tcl=-0.1, labels=FALSE)
        dev.off()
        #library(plyr)
        #library(randomforest)
        #library(pracma)
      }
    }
  }
}
