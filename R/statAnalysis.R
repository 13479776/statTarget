#' statAnalysis provide the statistical analysis for metabolomics data or others.
#' @param file a file with  the expression information. 
#' @param Frule The cut-off value for missing value filter function.
#' @param imputeM The parameter for imputation method.(i.e., nearest neighbor averaging, "KNN"; minimum values for imputed variables, "min", median values for imputed variables (Group dependent) "median"). 
#' @param glog is a data index for log transformation, with the default value TRUE.
#' @param test.multi must be an index for statistic analysis, with the default value TRUE.
#' @param scaling is index of scaling method index for statistic analysis (PCA or PLS-DA). 'pareto', 'Pareto', 'p' or 'P' can be used for specifying the Pareto scaling. 'auto', 'Auto', 'auto', 'a' or 'A' can be used for specifying the Auto scaling (or unit variance scaling). 'vast', 'Vast', 'v' or 'V' can be used for specifying the vast scaling. 'range', 'Range', 'r' or 'R' can be used for specifying the Range scaling.
#' @param nvarRF shows the number of variables in Gini plot of Randomforest model (=< 100). 
#' @param silt shows the number of permutation times
#' @param pcax shows principal components in PCA model for the x-axis.
#' @param pcay shows principal components in PCA model for the y-axis.
#' @export
statAnalysis <- function (file, Frule = 0.8,imputeM = "KNN", glog = TRUE, test.multi=TRUE, nvarRF =10, scaling = "Pareto",silt = 500, pcax = 1, pcay = 2) {
  dirout.uni = paste(getwd(), "/statTarget/", sep = "")
  dir.create(dirout.uni)
  dirout.uni = paste(getwd(), "/statTarget/statAnalysis/", sep = "")
  dir.create(dirout.uni)
  dat <- read.csv(file,header=TRUE)
  cat(date(), "\nstatTarget: statistical analysis start... \nEvaluation of missing value...")
  dat <- as.matrix(dat)
  dat[dat<=0] <- NA
  message(date(), "The number of NA value in Data Profile: ",
          sum(is.na(dat) | as.matrix(dat) <= 0))
  imdat <- dat
  #############Filter miss value###################
  
  FilterMV = function(m,degree) {
    dx <- c() 
    for(i in 1:ncol(m)){
      freq <- as.vector(tapply(m[,i], degree, function(x){sum(is.na(x) | as.matrix(x) <= 0)/length(x)}))
      if(sum(freq > Frule) > 0) dx <- c(dx , i)
    } 
    if(length(dx) >0) m <- m[,-dx]
    return(m)
  }
  classF <- as.factor(imdat[,2])
  #classF = addNA(classF)
  imdatF = FilterMV(imdat,classF)
  Frule_warning= paste("The number of vaiables including", Frule*100, "% of missing value :",sep = " ")
  message(Frule_warning," ", dim(imdat)[2]-dim(imdatF)[2])
  #imsamFP = imdatF
  #degree <- as.factor(dat[,2])
  #rule80 <- rule80(dat,degree)
  
  ##############impute missing value#################
  cat(date(), "\nImputation start...\n")
  if(imputeM == "KNN"){
    #require(impute)
    mvd <- impute::impute.knn(imdatF[,3:ncol(imdatF)])
    inputedData <- mvd$data
  }else if(imputeM == "min"){
    inputedData <- apply(imdatF[,3:ncol(imdatF)],2,function(y){
      y[is.na(y) | y<=0] <- min(y[y>0],na.rm = TRUE)
      y})    
    #inputedData <- t(inputedData)
  }else if(imputeM == "median"){
    missvalue <- function(x,group) {
      x[is.na(x) == TRUE ] <- 0
      group = as.factor(as.numeric(group))
      for (i in 1:dim(x)[1]){  
        for(j in 3:dim(x)[2]){
          if(x[i,j] == 0 | sum(is.na(x[i,j])) >= 1){
            #x[i,j][is.na(x[i,j]) == TRUE ] <- 0
            x[i,j] <- tapply(as.numeric(x[,j]),group,median)[group[i]]
          }
        }
      } 
      return(x)
    }
    inputedData = missvalue(imdatF,classF)
    inputedData = inputedData[,-c(1,2)]
  }
  message(date(), "The number of NA value in Data Profile after the initial imputation: ",
          sum(is.na(inputedData) | as.matrix(inputedData) <= 0))
  
  if(sum(is.na(inputedData) | as.matrix(inputedData) <= 0) > 0)
  {  
    inputedData[inputedData<=0] <- NA
    mvd2 <- impute::impute.knn(inputedData[,1:ncol(inputedData)])
    inputedData <- mvd2$data
    message(date(),"The number of NA value in Data Profile after the second imputation (KNN): ",
            sum(is.na(inputedData) | as.matrix(inputedData) <= 0))
  }
  
  cat("\nImputation Finished!\n")
  #msva <- missvalue(rule80,degree)
  dirout.uni = paste(getwd(), "/statTarget/statAnalysis/PreTable/", sep = "")
  dir.create(dirout.uni)
  prefile = paste(getwd(), "/statTarget/statAnalysis/PreTable/imputation.csv", sep = "")
  
  write.csv(cbind(imdatF[,1:2],inputedData), prefile, row.names = FALSE)
  bStatX(prefile)
  message(date(), "\\Statistic Summary Finished")
  if (glog) {
    #glog trans
    x <- read.csv(prefile, sep = ",", header = TRUE)
    GloggedSmpd<-glog(x[,3:ncol(x)],2)
    #GloggedSmpd[1:5,1:4]
    
    sdv <- apply(GloggedSmpd,1,sd)
    meanI <- apply(GloggedSmpd,1,mean)
    logvarI <- data.frame(meanI,sdv)
    log_rankI <- logvarI[order(logvarI[,1],decreasing=F),]
    
    sdvf <- apply(x[,3:ncol(x)],1,sd)
    meanII <- apply(x[,3:ncol(x)],1,mean)
    logvarII <- data.frame(meanII,sdvf)
    log_rankII <- logvarII[order(logvarII[,1],decreasing=F),]
    logfile = paste(getwd(), "/statTarget/statAnalysis/PreTable/Table_imputation_glog.csv", sep = "")
    write.csv(cbind(x[,1:2],GloggedSmpd), logfile, row.names = FALSE)
    #pdf(paste(getwd(), "/statTarget/hist_plot_glog.pdf"))
    pdf("./statTarget/statAnalysis/PreTable/hist_plot_glog.pdf")
    par(mfrow=c(1,2))
    plot(1:dim(log_rankI)[1],log_rankI[,2],pch= 21,bg="red",col=rgb(0,0,0,100,maxColorValue=255),xlab="rank of mean intensity", ylab="standard deviation")
    plot(1:dim(log_rankII)[1],log_rankII[,2],pch= 21,bg="green",col=rgb(0,0,0,100,maxColorValue=255),xlab="rank of mean intensity",ylab="standard deviation")
    dev.off()
    message(date(), "\\Preglog Finished!")
  } else {
    message(date(), "\\Preglog NONE!")
  }
  
  if (test.multi) {
    if (glog) {
      setwd("./statTarget/statAnalysis/")
      #par(mfrow=c(1,1))
      message(date(), "\\PCA-PLSDA start...")
      logf <-read.csv(logfile,header = TRUE)
      if(nvarRF > ncol(logf)-2)
      {
        stop("Do not set the value of nvarRF higher than the number of varibles")
      }
      explore.data.stat(logfile,scaling,normalize = TRUE)
      Plot.pca.score.stat(pcax,pcay,scaling)
      Plot.pca.loading(pcax,pcay,scaling)
      outlier.stat(pcax,pcay,scaling)
      message(date(), "\\Variable List")
      chose.driver(scaling)
      #message(date(),"\\PCA-PLSDA Start...")
      plsda.stat(scaling,silt)
      Plot.plsda.stat(1,2,scaling)#In par(new = T) : calling par(new=TRUE) with no plot
      log = sT.univariate(logfile,normalize = TRUE,nvarRF)
      oplsda(scaling)
      #dev.off()
      par("mfrow")
      message(date(), "\\Multiglog Finished!")
    }else {
      setwd("./statTarget/statAnalysis/")
      #par(mfrow=c(1,1))
      message(date(), "\\PCA-PLSDA start...")
      logFF <-read.csv(prefile,header = TRUE)
      if(nvarRF > ncol(logFF)-2)
      {
        stop("Do not set the value of nvarRF higher than the number of varibles")
      }
      
      explore.data.stat(prefile,scaling,normalize = TRUE)
      Plot.pca.score.stat(pcax,pcay,scaling)
      Plot.pca.loading(pcax,pcay,scaling)
      outlier.stat(pcax,pcay,scaling)
      message(date(), "\\Variable List")
      chose.driver(scaling)
      #message(date(), "\\PCA-PLSDA Start...")
      plsda.stat(scaling,silt)
      Plot.plsda.stat(1,2,scaling)
      log = sT.univariate(prefile,normalize = TRUE,nvarRF)
      oplsda(scaling)
      #dev.off()
      message(date(), "\\Multiglog FREE!")
    }
    message(date(), "\\Ttest.multi done!")
    cat(date(),"\\Statistical Analysis Finished!")
    tmpfile = paste(getwd(), "/tmp/",sep="")
    unlink(tmpfile, recursive=TRUE)
    tmpfiles = paste(getwd(), "/Groups/",sep="")
    unlink(tmpfiles, recursive=TRUE)
  }
}
