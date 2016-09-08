#' @name statAnalysis
#' @title statAnalysis for GUI
#' @description statAnalysis provide the statistical analysis for metabolomics
#'  data or others.
#' @param file The file with  the expression information. 
#' @param Frule The cut-off value for missing value filter function.
#' @param imputeM The parameter for imputation method.(i.e., nearest neighbor 
#' averaging, "KNN"; minimum values for imputed variables, "min", median values
#'  for imputed variables (Group dependent) "median"). 
#' @param glog  Generalised logarithm (glog) transformation, with the default 
#' value TRUE.
#' @param test.multi Multiple statistical analysis, with the default value TRUE.
#' @param FDR The false discovery rate for conceptualizing the rate of type I 
#' errors in null hypothesis testing when conducting multiple comparisons.
#' @param scaling Scaling method before statistic analysis (PCA or PLS-DA). 
#' 'pareto', 'Pareto', 'p' or 'P' can be used for specifying the Pareto scaling.
#'  'auto', 'Auto', 'auto', 'a' or 'A' can be used for specifying the Auto 
#'  scaling (or unit variance scaling). 'vast', 'Vast', 'v' or 'V' can be used
#'   for specifying the vast scaling. 'range', 'Range', 'r' or 'R' can be used 
#'   for specifying the Range scaling.
#' @param nvarRF The number of variables in Gini plot of Randomforest model
#'  (=< 100). 
#' @param silt The number of permutation times for PLS-DA model
#' @param pcax Principal components in PCA model for the x-axis.
#' @param pcay Principal components in PCA model for the y-axis.
#' @param Labels Name labels for score plot of multiple statistical analysis
#' @param upper.lim The up-regulated metabolites using Fold Changes cut off 
#' values in the Volcano plot.
#' @param lower.lim The down-regulated metabolites using Fold Changes cut off
#'  values in the Volcano plot.
#' @param sig.lim The significance level for metabolites in the Volcano plot.
#' @return A object of statAnalysis
#' @examples 
#' datpath <- system.file("extdata",package = "statTarget")
#' file <- paste(datpath,"data_example.csv", sep="/")
#' statAnalysis(file,nvarRF =5)
#' @export
#' @keywords PCA PLSDA P-value 
statAnalysis <- function (file, Frule = 0.8,imputeM = "KNN", glog = TRUE, 
                          test.multi=TRUE, FDR = TRUE, nvarRF =10, 
                          scaling = "Pareto",silt = 500, pcax = 1, pcay = 2,
                          Labels = TRUE, upper.lim = 1.5, lower.lim = 0.5, 
                          sig.lim = 0.05) {
  dirout.uni = paste(getwd(), "/statTarget/", sep = "")
  dir.create(dirout.uni)
  dirout.uni = paste(getwd(), "/statTarget/statAnalysis/", sep = "")
  dir.create(dirout.uni)
  dat <- read.csv(file,header=TRUE)
  cat("\n\nstatTarget: statistical analysis start... Time: ", date(),  
      "\n\nStep 1: Evaluation of missing value...")
  #dat <- as.matrix(dat)
  #dat[dat==0] <- NA
  message( "\nThe number of NA value in Data Profile: ",
          sum(is.na(dat) | as.matrix(dat) == 0))
  imdat <- dat
  #############Filter miss value###################
  
  FilterMV = function(m,degree) {
    dx <- c() 
    for(i in 1:ncol(m)){
      freq <- as.vector(tapply(m[,i], degree, function(x){
        sum(is.na(x) | as.matrix(x) == 0)/length(x)}))
      if(sum(freq > Frule) > 0) dx <- c(dx , i)
    } 
    if(length(dx) >0) m <- m[,-dx]
    return(m)
  }
  classF <- as.factor(imdat[,2])
  #classF = addNA(classF)
  imdatF = FilterMV(imdat,classF)
  Frule_warning= paste("\nThe number of vaiables including", Frule*100, 
                       "% of missing value :",sep = " ")
  message(Frule_warning," ", dim(imdat)[2]-dim(imdatF)[2])
  #imsamFP = imdatF
  #degree <- as.factor(dat[,2])
  #rule80 <- rule80(dat,degree)
  imdatF <- as.matrix(imdatF)
  ##############impute missing value#################
  cat("\nStep 2: Imputation start... Time: ", date())
  if(imputeM == "KNN"){
    #require(impute)
    mvd <- impute::impute.knn(imdatF[,3:ncol(imdatF)],rowmax = 0.99, 
                              colmax = 0.99, maxp = 15000)
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
  message( "\nThe number of NA value in Data Profile after the initial", 
           "imputation: ",
          sum(is.na(inputedData)| as.matrix(inputedData) == 0))
  
  if(sum(is.na(inputedData) | as.matrix(inputedData) == 0) > 0)
  {  
    inputedData[inputedData==0] <- NA
    mvd2 <- impute::impute.knn(inputedData[,1:ncol(inputedData)],
                               rowmax = 0.99, colmax = 0.99, maxp = 15000)
    inputedData <- mvd2$data
    message( "\nThe number of NA value in Data Profile after", 
             "the second imputation (KNN): ",
            sum(is.na(inputedData) | as.matrix(inputedData) == 0))
  }
  
  message("\nImputation Finished!")
  
  ##.....................................
  
  ##    Transform the Factor     
  
  ##.....................................
  
  TraceFc <- function(x){
    xF <- factor(x[,2])
    for(i in 1:length(levels(xF))) {
      levels(xF)[i] <- i
    }
    x[,2] <- xF
    return(x)
  }
  
  imdatF2 <- TraceFc(imdatF)
  
  dirout.uni = paste(getwd(), "/statTarget/statAnalysis/PreTable/", sep = "")
  dir.create(dirout.uni)
  
  mach <- cbind(data.frame(dat[,2]),data.frame(imdatF2[,2]))
  colnames(mach) <- c("Class","Number")
  machfile = paste(getwd(), "/statTarget/statAnalysis/PreTable/slink.csv", 
                   sep = "")
  write.csv(mach, machfile, row.names = FALSE)
  
  prefile = paste(getwd(), "/statTarget/statAnalysis/PreTable/imputation.csv",
                  sep = "")
  write.csv(cbind(imdatF2[,1:2],inputedData), prefile, row.names = FALSE)
  
  cat("\nStep 3: Statistic Summary Start... Time: ", date())
  bStatX(prefile)
  
  if (glog) {
    #glog trans
    x <- read.csv(prefile, sep = ",", header = TRUE)
    GloggedSmpd<-glog(x[,3:ncol(x)],2)
    #GloggedSmpd[1:5,1:4]
    
    sdv <- apply(GloggedSmpd,1,sd)
    meanI <- apply(GloggedSmpd,1,mean)
    logvarI <- data.frame(meanI,sdv)
    log_rankI <- logvarI[order(logvarI[,1],decreasing=FALSE),]
    
    sdvf <- apply(x[,3:ncol(x)],1,sd)
    meanII <- apply(x[,3:ncol(x)],1,mean)
    logvarII <- data.frame(meanII,sdvf)
    log_rankII <- logvarII[order(logvarII[,1],decreasing=FALSE),]
    logfile = paste(getwd(), 
              "/statTarget/statAnalysis/PreTable/Table_imputation_glog.csv", 
              sep = "")
    write.csv(cbind(x[,1:2],GloggedSmpd), logfile, row.names = FALSE)
    #pdf(paste(getwd(), "/statTarget/hist_plot_glog.pdf"))
    pdf("./statTarget/statAnalysis/PreTable/hist_plot_glog.pdf")
    par(mfrow=c(1,2))
    plot(1:dim(log_rankI)[1],log_rankI[,2],pch= 21,bg="red",
         col=rgb(0,0,0,100,maxColorValue=255),xlab="rank of mean intensity", 
         ylab="standard deviation")
    plot(1:dim(log_rankII)[1],log_rankII[,2],pch= 21,
         bg="green",col=rgb(0,0,0,100,maxColorValue=255),
         xlab="rank of mean intensity",ylab="standard deviation")
    dev.off()
    #message( "\nPreglog Finished!")
  } else {
    message( "\nPre-glog NONE!")
  }
  
  ##.....................................
  
  ##    Multiple statistical analysis     
  
  ##.....................................
  
  if (test.multi) {
    if (glog) {
      setwd("./statTarget/statAnalysis/")
      #par(mfrow=c(1,1))
      cat("\n\nStep 4: Glog PCA-PLSDA start... Time: ", date())
      logf <-read.csv(logfile,header = TRUE)
      if(nvarRF > ncol(logf)-2)
      {
        stop(
          "\nDo not set the value of nvarRF higher than the number of varibles")
      }
      explore_data_stat(logfile,scaling,normalize = TRUE)
      Plot_pca_score_stat(pcax,pcay,scaling,Labels)
      Plot_pca_loading(pcax,pcay,scaling)
      outlier_stat(pcax,pcay,scaling)
      #message( "\\Variable List")
      #chose.driver(scaling)
      #message(date(),"\\PCA-PLSDA Start...")
      plsda_stat(scaling,silt)
      Plot_plsda_stat(1,2,scaling,Labels)
      #In par(new = T) : calling par(new=TRUE) with no plot
      
      cat("\nStep 5: Univariate Test Start...! Time: ", date())
      log = sT_univariate(logfile,normalize = TRUE,
                          FDR=FDR, nvarRF = nvarRF,upper.lim = upper.lim, 
                          lower.lim = lower.lim, sig.lim = sig.lim)
      #pca.osc(scaling,Labels)
      #dev.off()
      #par("mfrow")
      #message("\\Multiglog Finished!", date())
    }else {
      setwd("./statTarget/statAnalysis/")
      #par(mfrow=c(1,1))
      cat("\n\nStep 4: PCA-PLSDA start... Time: ", date())
      logFF <-read.csv(prefile,header = TRUE)
      if(nvarRF > ncol(logFF)-2)
      {
        stop(
          "Do not set the value of nvarRF higher than the number of varibles")
      }
      
      explore_data_stat(prefile,scaling,normalize = TRUE)
      Plot_pca_score_stat(pcax,pcay,scaling,Labels)
      Plot_pca_loading(pcax,pcay,scaling)
      outlier_stat(pcax,pcay,scaling)
      #message(date(), "\\Variable List")
      #chose.driver(scaling)
      #message(date(), "\\PCA-PLSDA Start...")
      plsda_stat(scaling,silt)
      Plot_plsda_stat(1,2,scaling, Labels)
      #pca.osc(scaling,Labels)
      cat("\nStep 5: Univariate Test Start...! Time: ", date())
      log = sT_univariate(prefile,normalize = TRUE,FDR=FDR, 
                          nvarRF = nvarRF,upper.lim = upper.lim, 
                          lower.lim = lower.lim, sig.lim = sig.lim)
      #dev.off()
      #message(date(), "\\Multiglog FREE!")
    }
    #message(, "\\Ttest.multi done!",date())
    message("\nStatistical Analysis Finished! Time: ",date())
    tmpfile = paste(getwd(), "/tmp/",sep="")
    unlink(tmpfile, recursive=TRUE)
    tmpfiles = paste(getwd(), "/Groups/",sep="")
    unlink(tmpfiles, recursive=TRUE)
  }
}
