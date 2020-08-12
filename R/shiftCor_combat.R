#' @name shiftCor_dQC
#' @title QC-free based signal correction
#' @description  shiftCor_dQC provides the QC-free based signal correction for 
#' large scale mass spectrometry-based omics data.
#' @param samPeno The file with the meta information including the sample name,
#'  batches, class and order (denoting other covariates besides batch). 
#' @param samFile The file with the expression information. 
#' @param Frule Modified n precent rule function. A variable will be kept if it has a non-zero value
#' for at least n precent of samples in any one group. 
#' (Default: 0.8)  
#' @param MLmethod 'ComBat' allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects.The function was revised accroding to 'sva' package (version = "3.8").
#' @param imputeM The parameter for imputation method i.e., nearest neighbor 
#' averaging, 'KNN'; minimum values, 'min'; Half of minimum values, 'minHalf'; 
#' median values, 'median'.
#' @param par.prior TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots.
#' @param mod.covariates TRUE indicates model matrix for outcome of interest and other covariates besides batch (Column 'order' denotes covariates the  in samPeno file).
#' @param batch.Num (Optional) NULL If given, will use the selected batch as a reference for batch adjustment.
#' @return the shiftCor files. See the details at https://stattarget.github.io
#' @examples 
#' datpath <- system.file('extdata',package = 'statTarget')
#' samPeno <- paste(datpath,'MTBLS79_dQC_sampleList.csv', sep='/')
#' samFile <- paste(datpath,'MTBLS79.csv', sep='/')
#' shiftCor_dQC(samPeno,samFile, Frule = 0.8, MLmethod = "Combat",mod.covariates = FALSE)
#' shiftCor_dQC(samPeno,samFile, Frule = 0.8, MLmethod = "Combat",mod.covariates = TRUE,batch.Num = 1)
#' @keywords Quality Controls,Correction
#' @export 
shiftCor_dQC <- function(samPeno, samFile, Frule = 0.8, imputeM = "KNN",
                         MLmethod = "Combat",par.prior = TRUE,prior.plots = FALSE,
                         mod.covariates =FALSE, batch.Num = NULL) {
    cat("\n")
    
    cat("statTarget: Signal Correction Start... Time:", date(), "\n\n")
    cat("* Step 1: Data File Checking Start..., Time: ", date(), "\n")
    
    cat("\n", "Data Link", "\n")
    cat(" metaFile:", samPeno, "\n")
    cat(" profileFile:", samFile, "\n")
    
    samPeno <- read.csv(samPeno, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    samPeno <- as.data.frame(samPeno)
    #samPeno <- plyr::arrange(samPeno, order)
    samFile <- read.csv(samFile, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
    samFile <- t(samFile)
    colnames(samFile) <- samFile[1, ]
    samFile <- as.data.frame(samFile[-1, ])
    rownames(samFile) <- samFile$name
    
    ############## Checking the input file#############
    
    
    ## .............................................
    
    ## Data Matching
    
    ## .............................................
    
    cat("\n", " ", dim(samPeno)[1], " Meta Samples vs ", dim(samFile)[1], " Profile samples", sep = "")
    cat("\n", "The Meta samples list (*NA, missing data from the Profile File)")
    cat("\n\n")
    mcdat <- samFile[, 1][match(samPeno[, 1], samFile[, 1])]
    print(as.vector(mcdat))
    
    ## .............................................
    if (any(is.na(mcdat))) {
        stop("\n", "Missing data from the Profile File! Check your data please!!")
    }
    
    if (dim(samFile)[1] - dim(samPeno)[1] > 0) {
        cat("\n", "Warning: The sample size in Profile File is larger than Meta File! ")
        cat("\n")
    } else if (dim(samFile)[1] - dim(samPeno)[1] < 0) {
        stop("\n", "The sample size in Profile File should be no less than Meta File!", " Check your data please!!")
    }
    
    samFP <- samFile[samPeno$sample, ]
    
    if (sum(is.na(samPeno$class)) <= 0) {
        warning("\n", "There were no QC samples in your data!")
    }
    
    samPeno_stat <- samPeno
    rownames(samPeno_stat) <- samPeno[, 1]
    #qc_seq <- rownames(samPeno_stat)
    #qc_seq_tmp <- grep("QC", qc_seq)
    
    
    
    #sam_seq_tmp <- samPeno_stat[-c(qc_seq_tmp), ]
    
    if (sum(is.na(samPeno_stat$batch)) > 0) {
        stop("\n", "There were missing values in batch! Check your data please!\n")
    }
    if (sum(is.na(samPeno_stat$class)) > 0) {
        stop("\n", "There were missing values (Unclassified data) in sample class! 
         Check your data please!\n")
    }
    
    
    #samPeno_stat$class[is.na(samPeno_stat$class)] <- "QC"
    data_stat = aggregate(samPeno_stat$class, by = list(Category = samPeno_stat$class), FUN = length)
    batch_stat = aggregate(samPeno_stat$class, by = list(Category = samPeno_stat$batch), FUN = length)
    colnames(data_stat) <- c("Class", "No.")
    colnames(batch_stat) <- c("Batch", "No.")
    cat("\n", "Meta-information:", "\n\n")
    print(data_stat)
    print(batch_stat)
    num_sam <- dim(samFile[, 2:ncol(samFile)])
    num_sam <- data.frame(num_sam)
    colnames(num_sam) <- c("no.")
    rownames(num_sam) <- c("Samples", "Metabolites")
    cat("\n", "Metabolic profile information:", "\n\n")
    print(num_sam)
    
    
    ################# data NA########
    samFP_temp <- as.matrix(samFP[, 2:ncol(samFP)])
    
    if (length(samFP_temp[samFP_temp < 0L]) > 0) {
        samFP_temp[samFP_temp < 0L] <- 0L
    }
    samFP_temp[samFP_temp == 0L] <- NA
    
    cat("\n")
    cat("* Step 2: Evaluation of Missing Value...", "\n")
    cat("\n", "The number of missing value before signal correction: ", sum(is.na(samFP_temp)))
    imsamFP <- samFP_temp
    ############# Filter miss value###################
    
    FilterMV = function(m, degree) {
        dx <- c()
        for (i in 1:ncol(m)) {
            freq <- as.vector(tapply(m[, i], degree, function(x) {
                sum(is.na(x))/length(x)
            }))
            if (sum(freq > 1 - Frule) > 0) 
                dx <- c(dx, i)
        }
        if (length(dx) > 0) 
            m <- m[, -dx]
        return(m)
    }
    classF <- as.factor(samPeno$class)
    #classF = addNA(classF)
    imsamFPF = FilterMV(imsamFP, classF)
    Frule_warning = paste("The number of filtered variables using the modified ", Frule * 100, "% rule :", 
        sep = " ")
    cat("\n", Frule_warning, " ", dim(imsamFP)[2] - dim(imsamFPF)[2], "\n")
    imsamFP = as.matrix(imsamFPF)
    ############## impute missing value#################
    cat("\n")
    cat("* Step 3: Imputation start...", "\n")
    
    
    if (imputeM == "KNN") {
        # require(impute)
        cat("\n", "The imputation method was set at 'KNN'")
        mvd <- impute::impute.knn(imsamFP, rowmax = 0.99, colmax = 0.99, maxp = 15000)
        inputedData <- mvd$data
    } else if (imputeM == "min") {
        cat("\n", "The imputation method was set at 'min'")
        
        minValue <- function(x, group) {
            group = as.factor(as.numeric(group))
            for (i in 1:dim(x)[1]) {
                for (j in 1:dim(x)[2]) {
                  if (is.na(x[i, j]) == TRUE) {
                    x[i, j] <- tapply(as.numeric(x[, j]), group, min, na.rm = TRUE)[group[i]]
                  }
                }
            }
            return(x)
        }
        inputedData = minValue(imsamFP, classF)
        # inputedData = inputedData[,-1]
        
    } else if (imputeM == "minHalf") {
        cat("\n", "The imputation method was set at 'minHalf'")
        
        minHalfValue <- function(x, group) {
            group = as.factor(as.numeric(group))
            for (i in 1:dim(x)[1]) {
                for (j in 1:dim(x)[2]) {
                  if (is.na(x[i, j]) == TRUE) {
                    x[i, j] <- tapply(as.numeric(x[, j]), group, min, na.rm = TRUE)[group[i]]/2
                  }
                }
            }
            return(x)
        }
        inputedData = minHalfValue(imsamFP, classF)
        # inputedData = inputedData[,-1]
        
    } else if (imputeM == "median") {
        
        medianvalue <- function(x, group) {
            group = as.factor(as.numeric(group))
            for (i in 1:dim(x)[1]) {
                for (j in 1:dim(x)[2]) {
                  if (is.na(x[i, j]) == TRUE) {
                    x[i, j] <- tapply(as.numeric(x[, j]), group, median, na.rm = TRUE)[group[i]]
                  }
                }
            }
            return(x)
        }
        
        cat("\n", "The imputation method was set at 'median'")
        
        inputedData = medianvalue(imsamFP, classF)
        # inputedData = inputedData[,-1]
    }
    cat("\n", "The number of missing value after imputation: ", sum(is.na(inputedData)))
    
    
    if (sum(is.na(inputedData)) > 0) {
        inputedData = minHalfValue(inputedData, classF)
        cat("\n", "The number of missing value after the second imputation: ", sum(is.na(inputedData)))
    }
    
    
    cat("\n", "Imputation Finished!", "\n", "\n")
    
    cat("* Step 4: Signal Correction Start... Time: ", date(), "\n")
    
    dat <- as.matrix(t(inputedData))
    #numX <- 1:dim(dat)[2]
    
    
    if (MLmethod == "Combat") {
        # require(randomForest)
        cat("\n", "The Signal Correction method was set at Combat", "\n")
        
       #suprr <- capture.output(
        if(mod.covariates) {
          batchinfo <- samPeno_stat$batch
          covar <- as.factor(samPeno_stat[,4])
          if(length(levels(covar)) - length(covar) == 0 ) stop("No. of covariates levels should be not equal to No. of samples")
          mod = model.matrix(~ covar, data=samPeno_stat)
          comba <- ComBat_qcFree(dat=dat, batch=batchinfo, mod=mod, 
                                 par.prior=par.prior, prior.plots=prior.plots,
                                 ref.batch = batch.Num)
        } else {
          batchinfo <- samPeno_stat$batch
        comba <- ComBat_qcFree(dat=dat, batch=batchinfo, mod=NULL, 
                              par.prior=par.prior, prior.plots=prior.plots,
                              ref.batch = batch.Num)
        }
       #)
       }
    
    combOutput <- as.data.frame(cbind(samPeno_stat$sample,samPeno_stat$class,t(comba)),
                               stringsAsFactors =FALSE)
    colnames(combOutput)[1] <- "sample";colnames(combOutput)[2] <- "class";
    
    combraw <- as.data.frame(cbind(samPeno_stat$sample,samPeno_stat$class,t(dat)),
                                stringsAsFactors =FALSE)
    colnames(combraw)[1] <- "sample";colnames(combraw)[2] <- "class";
    
    
    ############### dataCheck

    dirout.uni = paste(getwd(), "/statTarget/", sep = "")
    dirsc.ID = getwd()
    dir.create(dirout.uni)
    dirout.w = paste(getwd(), "/statTarget/shiftCor", sep = "")
    dir.create(dirout.w)
    dirout.Bs = paste(getwd(), "/statTarget/shiftCor/Before_shiftCor", sep = "")
    dir.create(dirout.Bs)
    dirout.As = paste(getwd(), "/statTarget/shiftCor/After_shiftCor", sep = "")
    dir.create(dirout.As)

    
    combafter = paste("shift_cor_Combat", ".csv", sep = "")
    write.csv(combOutput, paste(dirout.As, combafter, sep = "/"), row.names = FALSE)
    
    combabefore = paste("shift_raw_Combat", ".csv", sep = "")
    write.csv(combraw, paste(dirout.Bs, combabefore, sep = "/"), row.names = FALSE)
    
    
    cat("", "Output Link:", getwd(), "\n")
    cat("\n", "Correction Finished! Time: ", date(), "\n")
    cat("\n", "####################################", "\n")
    cat(" # Software Version: statTarget 2.0 #", "\n")
    cat(" ####################################", "\n")
    
    ################## Loess Plot########################
    setwd(dirsc.ID)
    
    # parameter output
    scPam1 <- c("Frule", "MLmethod", "imputeM", "par.prior","mod.covariates","batch.Num")
    scPam2 <- c(Frule, MLmethod, imputeM, par.prior,mod.covariates, paste("No.", batch.Num,sep = ""))
    scpam <- data.frame(scPam1, scPam2)
    colnames(scpam) <- c("parameter", "value")
    par_sh = paste("statTarget/ParameterShiftCor", ".log", sep = "")
    write.table(scpam, paste(getwd(), par_sh, sep = "/"), row.names = FALSE)
    
    tmpfilesc = paste(getwd(), "/tmp", sep = "")
    unlink(tmpfilesc, recursive = TRUE)
}




rowVars <- function(x, ...) 
{
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x, ...)), ...)/(n - 1))
}

aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- apply(!is.na(sdat),1,sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat,g.bar,n,d.old,t2)
    sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
    d.new <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  #cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}

int.eprior <- function(sdat,g.hat,d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]		
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x),length(g),n,byrow=T)
    resid2 <- (dat-g)^2
    sum2 <- resid2%*%j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star,sum(g*LH)/sum(LH))
    d.star <- c(d.star,sum(d*LH)/sum(LH))
    #if(i%%1000==0){cat(i,'\n')}
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  adjust	
} 

