#' @name shiftCor
#' @title shiftCor for GUI
#' @description  shiftCor provide the QC-RLS correction for 
#' large scale metabolomics.
#' @param samPeno The file with  the meta information including the sample name,
#'  batches, class and order. 
#' @param samFile The file with  the expression information. 
#' @param Frule The cut-off value for missing value filter function.
#' @param QCspan The smoothing parameter which controls the bias-variance 
#' tradeoff. if the QCspan is set at '0', the generalised cross-validation 
#' will be performed to avoid overfitting the observed data.
#' @param degree Lets you specify local constant regression (i.e., the 
#' Nadaraya-Watson estimator, degree=0), local linear regression (degree=1), 
#' or local polynomial fits (degree=2, the default).
#' @param imputeM The parameter for imputation method.(i.e., nearest neighbor 
#' averaging, "KNN"; minimum values for imputed variables, "min", median values 
#' for imputed variables (Group dependent) "median"). 
#' @return A objects of shiftCor
#' @examples 
#' datpath <- system.file("extdata",package = "statTarget")
#' samPeno <- paste(datpath,"MTBLS79_sampleList.csv", sep="/")
#' samFile <- paste(datpath,"MTBLS79.csv", sep="/")
#' shiftCor(samPeno,samFile)
#' @keywords QC-RLS correction
#' @export 
shiftCor <- function(samPeno,samFile,Frule = 0.8,QCspan = 0.75, 
                    degree = 2,imputeM = "KNN"){
  cat("\nData File Checking Start..., Time: ",date(),"\n")
  samPeno <- read.csv(samPeno, header=TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE)
  samPeno <- as.data.frame(samPeno)
  samFile <- read.csv(samFile,header=FALSE, check.names = FALSE,
                      stringsAsFactors = FALSE)
  samFile <- t(samFile)
  colnames(samFile) <- samFile[1,]
  samFile <- as.data.frame(samFile[-1,])
  rownames(samFile) <- samFile$name
  
  ############## Checking the input file#############
  
  
  ##.............................................
  
               ##Data Matching
  
  ##.............................................
  
  message("\n",dim(samPeno)[1]," Pheno Samples x",dim(samFile)[1], 
          " Profile samples",sep="")
  message("\nThe Pheno samples list (*NA, missing data from the Profile File)")
  mcdat <- samFile[,1][match(samPeno[,1],samFile[,1])]
  print(as.vector(mcdat))
  
  ##.............................................
  if(any(is.na(mcdat))) {
    stop("\nMissing data from the Profile File! Check your data please!!")
  } 
  
  if(dim(samFile)[1] - dim(samPeno)[1]>0){
    message(
      "\nWarning: The sample size in Profile File is larger than Pheno File! ")
  }else if(dim(samFile)[1] - dim(samPeno)[1]<0) {
    stop(
      "\nThe sample size in Profile File should be no less than Pheno File!
      Check your data please!!")
  }
  
  samFP <- samFile[samPeno$sample,]
  
  if(sum(is.na(samPeno$class)) <=0 ){
    stop("\nThere were not QC sample in your data!")
  }
  
  samPeno_stat <- samPeno
  rownames(samPeno_stat) <- samPeno[,1]
  qc_seq <- rownames(samPeno_stat)
  qc_seq_tmp <- grep("QC",qc_seq)
  
  
  
  sam_seq_tmp <- samPeno_stat[-c(qc_seq_tmp),]

  if(sum(is.na(samPeno_stat$batch)) > 0 ){
    stop("\nThere were NA values in batch! Check your data please!\n")
  }
  if(sum(is.na(sam_seq_tmp$class)) > 0 ){
    stop("\nThere were NA values (ungrouped data) in sample class! 
         Check your data please!\n")
  }
  
  
  samPeno_stat$class[is.na(samPeno_stat$class)] <- "QC"
  data_stat =  aggregate(samPeno_stat$class, 
                         by=list(Category=samPeno_stat$class), FUN=length)
  batch_stat = aggregate(samPeno_stat$class, 
                         by=list(Category=samPeno_stat$batch), FUN=length)
  colnames(data_stat) <-c("Class","No.")
  colnames(batch_stat) <-c("Batch","No.")
  message("\nPheno information:")
  print(data_stat)
  print(batch_stat)
  num_sam <- dim(samFile[,2:ncol(samFile)])
  num_sam <- data.frame(num_sam)
  colnames(num_sam) <- c("No.")
  rownames(num_sam) <- c("QC and samples","Metabolites")
  message("\nProfile information:")
  print(num_sam)
  
  
  #################LOESS data NA########
  samFP <- as.matrix(samFP)
  samFP[samFP==0] <- NA
  cat("\nstatTarget: shiftCor start...Time: ",date(), 
      "\n\nStep 1: Evaluation of missing value...")
  message("\nThe number of NA value in Data Profile before QC-RLSC: ",
          sum(is.na(samFP) | as.matrix(samFP) == 0))
  imsamFP <- samFP
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
  classF <- as.factor(samPeno$class)
  classF = addNA(classF)
  imsamFPF = FilterMV(imsamFP,classF)
  Frule_warning= paste("\nThe number of vaiables including", 
                       Frule*100, "% of missing value :",sep = " ")
  message(Frule_warning," ", dim(imsamFP)[2]-dim(imsamFPF)[2])
  imsamFP = as.matrix(imsamFPF)
  ##############impute missing value#################
  cat( "\nStep 2: Imputation start...\n")
  

  if(imputeM == "KNN"){
    #require(impute)
    mvd <- impute::impute.knn(imsamFP[,2:ncol(imsamFP)],
                              rowmax = 0.99, colmax = 0.99, maxp = 15000)
    inputedData <- mvd$data
  }else if(imputeM == "min"){
    inputedData <- apply(imsamFP[,2:ncol(imsamFP)],2,function(y){
      y[is.na(y) | y<=0] <- min(y[y>0],na.rm = TRUE)
      y})    
    #inputedData <- t(inputedData)
  }else if(imputeM == "median"){
    missvalue <- function(x,group) {
      x[is.na(x) == TRUE ] <- 0
      group = as.factor(as.numeric(group))
      for (i in 1:dim(x)[1]){  
        for(j in 2:dim(x)[2]){
          if(x[i,j] == 0 | sum(is.na(x[i,j])) >= 1){
            #x[i,j][is.na(x[i,j]) == TRUE ] <- 0
            x[i,j] <- tapply(as.numeric(x[,j]),group,median)[group[i]]
          }
        }
      } 
      return(x)
     }
    inputedData = missvalue(imsamFP,classF)
    inputedData = inputedData[,-1]
  }
  message(
    "\nThe number of NA value in Data Profile after the initial imputation: ",
          sum(is.na(inputedData) | as.matrix(inputedData) == 0))
  
  if(sum(is.na(inputedData) | as.matrix(inputedData) == 0) > 0)
  {  
    inputedData[inputedData==0] <- NA
    mvd2 <- impute::impute.knn(inputedData[,1:ncol(inputedData)],
                               rowmax = 0.99, colmax = 0.99, maxp = 15000)
    inputedData <- mvd2$data
    message(
 "\nThe number of NA value in Data Profile after the second imputation (KNN): ",
            sum(is.na(inputedData) | as.matrix(inputedData) == 0))
  }
  

  message("\nImputation Finished!")
          
  cat("\nStep 3: QC-RLSC Start... Time: ", date())
  
  dat <- as.matrix(t(inputedData))
  numX <- 1:dim(dat)[2]
  
  if(QCspan > 0){
    message("\nWarning: The QCspan was set at ", QCspan,"\n")
    
    loessFit=function(x,y,QCspan,degree){
      cn <- colnames(x)
      ####Check########
      st_QC<-grep("QC",cn[1])
      ed_QC<-grep("QC",cn[length(cn)])
      if(length(st_QC)==0)
      {
        stop(
          "\nWrong: the first sample must be QC sample; please check ......");
      }
      if(length(ed_QC)==0)
      {
        stop(
          "\nWrong: the sample at the end of sequence must be QC sample; 
          please check ......");
      }
      qcid <- grep("QC",cn)
      
      pb <- txtProgressBar(min = 1, max = dim(x)[1], style = 3)
      
      for(i in 1:dim(x)[1]){
        loe <- stats::loess(x[i,qcid]~qcid,span=QCspan,degree = degree)
        yf <- stats::predict(loe,y)
        x[i,] <- as.numeric(x[i,])/yf
        setTxtProgressBar(pb, i)
      }
      close(pb)
      loessDat = x
    }
    loessDat <- loessFit(x = dat,y=numX,QCspan = QCspan,degree = degree)
  }else if(QCspan <= 0){
    message(
      "\nWarning: The QCspan was set at '0'.\n", 
      "\nThe LOESS based generalised cross-validation was used to", 
      "avoid overfitting the observed data\n")
    autoFit <- function(xl,y){
      cn <- colnames(xl)
      ####Check########
      st_QC<-grep("QC",cn[1])
      ed_QC<-grep("QC",cn[length(cn)])
      if(length(st_QC)==0)
      {
        stop("Wrong: the first sample must be QC sample; please check ......");
      }
      if(length(ed_QC)==0)
      {
        stop(
          "Wrong: the sample at the end of sequence must be QC sample; 
          please check ......");
      }
      qcid <- grep("QC",cn)
      pb <- txtProgressBar(min = 1, max = dim(xl)[1], style = 3)
      
      for(i in 1:dim(xl)[1]){
        
        Sys.sleep(0.000001)
        
        loe1 <- loess(xl[i,qcid]~qcid)
        #loe2 <- loe1
        env <- environment() 
        sploe <- function(sp){
          loe2 <- get("loe1", envir=env)
          mod <- stats::update(loe2, span = sp)
          CVspan = loessGCV(mod)[["gcv"]]
        }
        sp <- c(seq(0.08,0.75,0.01))
        CVspan = as.matrix(lapply(sp,sploe))
        CVspan[!is.finite(as.numeric(CVspan))] <- NA
        minG <- data.frame(sp,CVspan)
        minspan <- minG[which.min(minG[,2]),1]
        minspan
        #sp <- c(seq(0.05,0.75,0.01))
        #CVspan <-c()
        #for(j in 1:length(sp)){
        #  mod <- stats::update(loe1, span = sp[j])
        #  CVspan[j] = loessGCV(mod)[["gcv"]]
        #}
        #minG <- as.matrix(data.frame(sp,CVspan))
        #minG[!is.finite(minG)] <- max(minG[,2],na.rm = TRUE)
        #minspan <- minG[which.min(minG[,2]),1]
        #minspan
        loeN <- stats::update(loe1, span = minspan)
        yf <- predict(loeN,y)
        xl[i,] <- as.numeric(xl[i,])/yf
        setTxtProgressBar(pb, i)
      }
      close(pb)
      loessDat = xl
    }
    loessDat <- autoFit(xl = dat,y = numX)
    
    
  }
 
###############
  loessDat1 <- apply(loessDat, 2,function(x) ifelse(x<0, 0, x))
  
  if(sum(is.na(loessDat) | as.matrix(loessDat) == 0) > 0)
  {  
    loessDat[loessDat==0] <- NA
    mvd2 <- impute::impute.knn(loessDat[,1:ncol(loessDat)], 
                               rowmax = 0.99, colmax = 0.99, maxp = 15000)
    loessDat <- mvd2$data
    message(
      "\nThe number of NA value in Data Profile after Loess Correction (KNN): ",
            sum(is.na(loessDat) | as.matrix(loessDat) == 0))
  }
  dirout.uni = paste(getwd(), "/statTarget/", sep = "")
  dir.create(dirout.uni)
  dirout.w = paste(getwd(), "/statTarget/shiftCor", sep="")
  dir.create(dirout.w)
  dirout.Bs = paste(getwd(), "/statTarget/shiftCor/Before_shiftCor", sep="")
  dir.create(dirout.Bs)
  dirout.As = paste(getwd(), "/statTarget/shiftCor/After_shiftCor", sep="")
  dir.create(dirout.As)
  ########### Out plot of each peak ############
  cat("\nHigh-resulution images output...")
          
  for(i in 1 :dim(dat)[1]){
    loplot(dat,loessDat,i)
  }
  ###############Raw output###########
  
  raw_temp <- cbind(samPeno,inputedData)
  nam_qc <- rownames(raw_temp)
  
  #####QC Cal################
  QC_temp_raw <- grep("QC",nam_qc)
  QC_temp_raw <- raw_temp[c(QC_temp_raw),]
  raw_temp_qc <- QC_temp_raw[,-c(3,4)]
  rownames(raw_temp_qc) <- NULL
  RSD30_CV=paste("shift_QC_raw",".csv", sep="")
  write.csv(raw_temp_qc, paste(dirout.Bs, RSD30_CV, sep="/"))
  cat("\n\nCalculation of CV distribution of raw peaks (QC)...\n\n")
  #raw_temp_qc = cbind(seq(1,dim(raw_temp_qc)[1],1),raw_temp_qc)
  Rsdist_QC_raw = RsdCal(raw_temp_qc,batch = TRUE,
                         DistPattern = TRUE, output = FALSE)
  
  #####Sample Cal################
  sam_temp_raw <- grep("QC",nam_qc)
  sam_temp_raw <- raw_temp[-c(sam_temp_raw),]
  raw_temp_sam <- sam_temp_raw[,-c(3,4)]
  rownames(raw_temp_sam) <- NULL
  RSD30_CV=paste("shift_sam_raw",".csv", sep="")
  write.csv(raw_temp_sam, paste(dirout.Bs, RSD30_CV, sep="/"))
  #Rsdist_sam_raw = RsdCal(raw_temp_sam,DistPattern = FALSE)
  #raw_temp_sam = cbind(seq(1,dim(raw_temp_sam)[1],1),raw_temp_sam)
  Rsdist_sam_raw = RsdCal(raw_temp_sam,batch = FALSE, 
                          DistPattern = FALSE, output = FALSE)
  
  
  ###############Loess output###########
  
  lo_temp <- cbind(samPeno,t(loessDat))
  nam_qc <- rownames(lo_temp)
  QC_temp <- grep("QC",nam_qc)
  QC_temp <- lo_temp[-c(QC_temp),]
  lo_temp_sam <- QC_temp[,-c(2,4)]
  rownames(lo_temp_sam) <- NULL
  RSD30_CV=paste("shift_sample_loess",".csv", sep="")
  write.csv(lo_temp_sam, paste(dirout.As, RSD30_CV, sep="/"))
  
  ############loess sample Cal#############
  #Rsdist_sam_cor = RsdCal(lo_temp_sam,DistPattern = FALSE)
  #lo_temp_sam = cbind(seq(1,dim(lo_temp_sam)[1],1),lo_temp_sam)
  Rsdist_sam_cor = RsdCal(lo_temp_sam,batch = FALSE, 
                          DistPattern = FALSE,output = FALSE)
  
  
  
  ############loess QC Cal#############
  QC_temp <- grep("QC",nam_qc)
  QC_temp <- lo_temp[c(QC_temp),]
  lo_temp_qc <- QC_temp[,-c(3,4)]
  rownames(lo_temp_qc) <- NULL
  RSD30_CV=paste("shift_QC_loess",".csv", sep="")
  write.csv(lo_temp_qc, paste(dirout.As, RSD30_CV, sep="/"))
  cat("\n\nCalculation of CV distribution of corrected peaks (QC)...\n\n")
  #Rsdist_QC_cor=RsdCal(lo_temp_qc,DistPattern = TRUE)
  #lo_temp_qc = cbind(seq(1,dim(lo_temp_qc)[1],1),lo_temp_qc)
  Rsdist_QC_cor=RsdCal(lo_temp_qc,batch = TRUE, 
                       DistPattern = TRUE,output = TRUE)
  
  #################
  RSDdist(Rsdist_sam_raw,Rsdist_sam_cor,Rsdist_QC_raw,Rsdist_QC_cor)
  lo_temp_all <- lo_temp[,-c(2,4)]
  rownames(lo_temp_all) <- NULL
  RSD30_CV=paste("shift_all_loess",".csv", sep="")
  write.csv(lo_temp_all, paste(dirout.As, RSD30_CV, sep="/"))
  cat("\n\nCorrection Finished! Time: ",date())
  ##################Loess Plot########################
}

