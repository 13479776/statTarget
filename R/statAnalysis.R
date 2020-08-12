#' @name statAnalysis
#' @title statAnalysis for statistical analysis for omics data or others.
#' @description statAnalysis provides the statistical analysis for metabolomics
#'  data or others.
#' @param file The file with the expression information. 
#' @param Frule Modified n precent rule function. A variable will be kept if it has a non-zero value
#' for at least n precent of samples in any one group. (Default: 0.8)  
#' @param imputeM The parameter for imputation method i.e., nearest neighbor 
#' averaging, 'KNN'; minimum values, 'min'; Half of minimum values, 'minHalf'; 
#' median values, 'median'. 
#' @param normM The parameter for normalization method (i.e median quotient 
#'  normalization, 'PQN'; integral normalization , 'SUM', and 'NONE').
#' @param glog  Generalised logarithm (glog) transformation, with the default 
#' value TRUE. The glog is a better behaved log transformation when some data 
#' values are zero or just near zero.
#' @param FDR The false discovery rate for conceptualizing the rate of type I 
#' errors in null hypothesis testing when conducting multiple comparisons.
#' @param scaling Scaling method before statistic analysis (PCA or PLS-DA). 
#' 'pareto', 'Pareto', 'p' or 'P' can be used for specifying the Pareto scaling.
#'  'auto', 'Auto', 'auto', 'a' or 'A' can be used for specifying the Auto 
#'  scaling (or unit variance scaling). 'vast', 'Vast', 'v' or 'V' can be used
#'   for specifying the vast scaling. 'range', 'Range', 'r' or 'R' can be used 
#'   for specifying the Range scaling.
#' @param ntree Number of trees to grow for randomForest model. This should not 
#' be set to too small a number, to ensure that every input row gets predicted 
#' at least a few times.
#' @param nvarRF The number of the variables with top importance in randomforest
#'  model
#' @param silt The number of permutation for PLS-DA model and variable importance of randomForest.
#' @param pcax Principal components in PCA model for the x-axis.
#' @param pcay Principal components in PCA model for the y-axis.
#' @param Labels Name labels for score plot of multiple statistical analysis
#' @param save.boxplot if TRUE, the box plot is performed
#' @param plot.volcano if TRUE, the volcano plot is performed
#' @param upper.lim The up-regulated metabolites using Fold Changes cut off 
#' values in the Volcano plot. Fold change values will be calculated 
#' before normalization step. 
#' @param lower.lim The down-regulated metabolites using Fold Changes cut off
#' values in the Volcano plot. Fold change values will be calculated before 
#' normalization step.
#' @param sig.lim The significance level for metabolites in the Pvalues file in 
#' the Volcano plot.
#' @return The statAnalsis output files. See the details at https://stattarget.github.io
#' @examples 
#' datpath <- system.file('extdata',package = 'statTarget')
#' file <- paste(datpath,'data_example.csv', sep='/')
#' statAnalysis(file,Frule = 0.8, normM = 'NONE', imputeM = 'KNN', glog = TRUE,scaling = 'Pareto')
#' @author Hemi Luan, hemi.luan@gmail.com
#' @export
#' @keywords PCA PLSDA P-value 
statAnalysis <- function(file, Frule = 0.8, normM = "NONE", imputeM = "KNN", glog = TRUE, FDR = TRUE, 
    ntree = 500, nvarRF = 5, scaling = "Pareto", plot.volcano = TRUE, save.boxplot =FALSE,silt = 20, pcax = 1, pcay = 2, Labels = TRUE, upper.lim = 2, 
    lower.lim = 0.5, sig.lim = 0.05) {
    dirout.uni = paste(getwd(), "/statTarget/", sep = "")
    dirsc.IDA = getwd()
    
    dir.create(dirout.uni)
    dirout.uni = paste(getwd(), "/statTarget/statAnalysis/", sep = "")
    dir.create(dirout.uni)
    dat <- read.csv(file, header = TRUE)
    cat("\n\nstatTarget: statistical analysis start... Time: ", date(), "\n\n")
    
    cat("* Step 1: Evaluation of missing value...", "\n")
    
    
    cat("\n", "Data Link", "\n")
    cat(" statFile:", file, "\n")
    
    imdat <- dat[, 3:ncol(dat)]
    
    if (length(imdat[imdat < 0L]) > 0) {
        imdat[imdat < 0L] <- 0L
    }
    imdat[imdat == 0L] <- NA
    
    imdat <- cbind(dat[, 1:2], imdat)
    
    cat("\n", "The number of missing value in Data Profile: ", sum(is.na(imdat)))
    
    ############# Filter miss value###################
    
    FilterMV = function(m, degree) {
        dx <- c()
        for (i in 1:ncol(m)) {
            freq <- as.vector(tapply(m[, i], degree, function(x) {
                sum(is.na(x) | as.matrix(x) == 0)/length(x)
            }))
            if (sum(freq > 1 - Frule) > 0) 
                dx <- c(dx, i)
        }
        if (length(dx) > 0) 
            m <- m[, -dx]
        return(m)
    }
    classF <- as.factor(imdat[, 2])
    imdatF = FilterMV(imdat, classF)
    Frule_warning = paste("The number of filtered variables using the modified ", Frule * 100, "% rule :", 
        sep = " ")
    cat("\n", Frule_warning, " ", dim(imdat)[2] - dim(imdatF)[2], "\n\n")
    imdatF <- data.frame(imdatF)
    
    dirout.uni = paste(getwd(), "/statTarget/statAnalysis/DataPretreatment/", sep = "")
    dir.create(dirout.uni)
    
    cat("* Step 2: Summary statistics start... Time: ", date(), "\n\n")
    
    
    cat("* Step 3: Missing value imputation start... Time: ", date(), "\n")
    
    ############## impute missing value#################
    cat("\n", "Imputation method was set at", imputeM)
    
    if (imputeM == "KNN") {
        # require(impute)
        mvd <- impute::impute.knn(as.matrix(imdatF[, 3:ncol(imdatF)]), rowmax = 0.99, colmax = 0.99, 
            maxp = 15000)
        inputedData <- mvd$data
    } else if (imputeM == "min") {
        
        minValue <- function(x, group) {
            group = as.factor(as.numeric(group))
            for (i in 1:dim(x)[1]) {
                for (j in 3:dim(x)[2]) {
                  if (is.na(x[i, j]) == TRUE) {
                    x[i, j] <- tapply(as.numeric(x[, j]), group, min, na.rm = TRUE)[group[i]]
                  }
                }
            }
            return(x)
        }
        inputedData = minValue(imdatF, classF)
        inputedData = inputedData[, -c(1, 2)]
        
    } else if (imputeM == "minHalf") {
        
        minHalfValue <- function(x, group) {
            group = as.factor(as.numeric(group))
            for (i in 1:dim(x)[1]) {
                for (j in 3:dim(x)[2]) {
                  if (is.na(x[i, j]) == TRUE) {
                    x[i, j] <- tapply(as.numeric(x[, j]), group, min, na.rm = TRUE)[group[i]]/2
                  }
                }
            }
            return(x)
        }
        inputedData = minHalfValue(imdatF, classF)
        inputedData = inputedData[, -c(1, 2)]
        
    } else if (imputeM == "median") {
        medianvalue <- function(x, group) {
            group = as.factor(as.numeric(group))
            for (i in 1:dim(x)[1]) {
                for (j in 3:dim(x)[2]) {
                  if (is.na(x[i, j]) == TRUE) {
                    x[i, j] <- tapply(as.numeric(x[, j]), group, median, na.rm = TRUE)[group[i]]
                  }
                }
            }
            return(x)
        }
        inputedData = medianvalue(imdatF, classF)
        inputedData = inputedData[, -c(1, 2)]
    }
    cat("\n", "The number of NA value after imputation:", sum(is.na(inputedData)))
    
    if (sum(is.na(inputedData)) > 0) {
        minHalfValue2 <- function(x, group) {
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
        inputedData = minHalfValue2(inputedData, classF)
        
        # mvd2 <- impute::impute.knn(as.matrix(inputedData[,1:ncol(inputedData)]), rowmax = 0.99, colmax =
        # 0.99, maxp = 15000) inputedData <- mvd2$data
        cat("\n", "The number of missing value after", "the second imputation: ", sum(is.na(inputedData)))
    }
    
    cat("\n", "Imputation Finished!", "\n")
    
    ## .....................................
    
    ## Transform the Factor
    
    ## .....................................
    
    TraceFc <- function(x) {
        xF <- factor(x[, 2])
        for (i in 1:length(levels(xF))) {
            levels(xF)[i] <- i
        }
        x[, 2] <- xF
        return(x)
    }
    
    imdatF2 <- TraceFc(imdatF)
    
    
    mach <- cbind(data.frame(dat[, 2]), data.frame(imdatF2[, 2]))
    colnames(mach) <- c("Class", "Number")
    machfile = paste(getwd(), "/statTarget/statAnalysis/DataPretreatment/slink.csv", sep = "")
    write.csv(mach, machfile, row.names = FALSE)
    
    prefile = paste(getwd(), "/statTarget/statAnalysis/DataPretreatment/data_imputation_", imputeM, 
        ".csv", sep = "")
    write.csv(cbind(imdatF[, 1:2], inputedData), prefile, row.names = FALSE)
    
    # Summary statistics
    imdatStat <- TraceFc(imdatF)
    bStatX(imdatStat)
    
    ############## Normalization #################
    cat("\n", "* Step 4: Normalization start... Time: ", date(), "\n")
    normTemp <- t(inputedData)
    normimdatF <- normTarget(normTemp, method = normM)
    
    cat("\n", "Normalization method was set at", normM, "\n")
    
    normF <- cbind(imdatF2[, 1:2], t(normimdatF))
    normF <- as.matrix(normF)
    # dataNorm output
    normfile = paste(getwd(), "/statTarget/statAnalysis/DataPretreatment/data_normalization_", normM, 
        ".csv", sep = "")
    write.csv(normF, normfile, row.names = FALSE)
    # dataPreprocess inputedData <- imdatF
    
    if (glog) {
        # glog trans
        x <- read.csv(normfile, header = TRUE)
        GloggedSmpd <- glog(x[, 3:ncol(x)], 2)
        sdv <- apply(GloggedSmpd, 2, sd)
        meanI <- apply(GloggedSmpd, 2, mean)
        logvarI <- data.frame(meanI, sdv)
        log_rankI <- logvarI[order(logvarI[, 1], decreasing = FALSE), ]
        
        sdvf <- apply(x[, 3:ncol(x)], 2, sd)
        meanII <- apply(x[, 3:ncol(x)], 2, mean)
        logvarII <- data.frame(meanII, sdvf)
        log_rankII <- logvarII[order(logvarII[, 1], decreasing = FALSE), ]
        
        # dataLog output
        logfile = paste(getwd(), "/statTarget/statAnalysis/DataPretreatment/data_glog.csv", sep = "")
        write.csv(cbind(x[, 1:2], GloggedSmpd), logfile, row.names = FALSE)
        pdf("./statTarget/statAnalysis/DataPretreatment/glogPlot.pdf")
        par(mfrow = c(1, 2))
        plot(1:dim(log_rankII)[1], log_rankII[, 2], pch = 21, bg = "green", col = rgb(0, 0, 0, 100, 
            maxColorValue = 255), xlab = "rank of mean intensity", ylab = "Standard deviation")
        plot(1:dim(log_rankI)[1], log_rankI[, 2], pch = 21, bg = "red", col = rgb(0, 0, 0, 100, maxColorValue = 255), 
            xlab = "rank of mean intensity", ylab = "Standard deviation (after glog-transformation)")
        graphics::title("Variance stabilization with glog-transformation", line = -3, outer = TRUE)
        dev.off()
    } else {
        cat("\n", "Warning: glog-transformation skipped!", "\n")
    }
    
    ## .....................................
    
    ## Multiple statistical analysis
    
    ## .....................................
    
    
    if (glog) {
        setwd("./statTarget/statAnalysis/")
        cat("\n")
        cat("* Step 5: Glog PCA-PLSDA start... Time: ", date(), "\n")
        
        cat("\n", "Scaling method was set at", scaling, "\n")
        explore_data_stat(logfile, scaling, normalize = FALSE)
        Plot_pca_score_stat(pcax, pcay, scaling, Labels)
        Plot_pca_loading(pcax, pcay, scaling)
        outlier_stat(pcax, pcay, scaling)
        plsda_stat(scaling, silt)
        Plot_plsda_stat(1, 2, scaling, Labels)
        cat("\n")
        
        cat("* Step 6: Glog Random Forest start... Time: ", date(), "\n")
        logf <- read.csv(logfile, header = TRUE)
        ST_rForest(logf, ntree = ntree, times = silt, nvarRF = nvarRF, Labels = Labels)
        
        cat("\n")
        cat("* Step 7: Univariate Test Start...! Time: ", date(), "\n")
        log = sT_univariate(logfile, FDR = FDR, upper.lim = upper.lim, lower.lim = lower.lim, sig.lim = sig.lim,plot.volcano = plot.volcano, save.boxplot = save.boxplot)
    } else {
        setwd("./statTarget/statAnalysis/")
        cat("\n", "Step 5: PCA-PLSDA start... Time: ", date(), "\n")
        cat("\n", "Scaling method was set at", scaling, "\n")
        explore_data_stat(normfile, scaling, normalize = FALSE)
        Plot_pca_score_stat(pcax, pcay, scaling, Labels)
        Plot_pca_loading(pcax, pcay, scaling)
        outlier_stat(pcax, pcay, scaling)
        plsda_stat(scaling, silt)
        Plot_plsda_stat(1, 2, scaling, Labels)
        cat("\n")
        cat("* Step 6: Random Forest start... Time: ", date(), "\n")
        logFF <- read.csv(normfile, header = TRUE)
        ST_rForest(logFF, ntree = ntree, times = silt, nvarRF = nvarRF, Labels = Labels)
        
        cat("\n")
        cat("* Step 7: Univariate Test Start...! Time: ", date(), "\n")
        log = sT_univariate(normfile, FDR = FDR, upper.lim = upper.lim, lower.lim = lower.lim, sig.lim = sig.lim, plot.volcano = plot.volcano, save.boxplot = save.boxplot)
    }
    cat("\n", "Output Link:", getwd(), "\n")
    cat("\n", "Statistical Analysis Finished! Time: ", date(), "\n")
    #cat("\n", "Correction Finished! Time: ", date(), "\n")
    cat("\n", "####################################", "\n")
    cat(" # Software Version: statTarget 2.0 #", "\n")
    cat(" ####################################", "\n")
    
    # Parameter set
    
    stPam1 <- c("Frule", "normM", "imputeM", "glog", "FDR", "ntree", "nvarRF", "scaling", "silt", 
        "pcax", "pcay", "Labels", "upper.lim", "lower.lim", "sig.lim")
    stPam2 <- c(Frule, normM, imputeM, glog, FDR, ntree, nvarRF, scaling, silt, pcax, pcay, Labels, 
        upper.lim, lower.lim, sig.lim)
    stpam <- data.frame(stPam1, stPam2)
    colnames(stpam) <- c("parameter", "value")
    par_st = paste("statTarget/ParameterStatAnalysis", ".log", sep = "")
    write.table(stpam, paste(dirsc.IDA, par_st, sep = "/"), row.names = FALSE)
    
    tmpfile = paste(getwd(), "/tmp/", sep = "")
    unlink(tmpfile, recursive = TRUE)
    tmpfiles = paste(getwd(), "/Groups/", sep = "")
    unlink(tmpfiles, recursive = TRUE)
    setwd(dirsc.IDA)
    
}
normTarget <- function(object, method = "PQN", use_percent = 100) {
    # Note features should be in rows, samples in columns Dieterle F,Ross A, Schlotterbeck G, Senn H.:
    # Probabilistic Quotient Normalization as Robust Method to Account for Diluition of Complex
    # Biological Mixtures. Application in 1H NMR Metabolomics.  Anal Chem 2006;78:4281-90.
    
    if (method == "PQN") {
        # SUM Normalization
        x.s <- matrix(colSums(object, na.rm = TRUE), nrow = 1)
        uni = matrix(rep(1, nrow(object)), ncol = 1)
        area.uni <- uni %*% x.s
        normSUM <- object/area.uni
        # PQN
        fectures <- abs(normSUM * 100)
        rem <- floor(nrow(fectures) - (nrow(fectures) * use_percent/100))
        # calculate variance
        varian <- apply(fectures, 1, stats::var)  #calculate variance for each feature across samples
        # remove features with highest variance
        sel <- order(varian)[1:(nrow(fectures) - rem)]
        fectures2 <- fectures[sel, ]  #data matrix that contains only features with low variance
        refer <- apply(fectures2, 1, median, na.rm = TRUE)
        quotient <- fectures2/refer
        quotient.median <- apply(quotient, 2, median, na.rm = TRUE)
        norm.fectures <- t(t(fectures2)/quotient.median)
        rownames(norm.fectures) <- rownames(fectures2)
        colnames(norm.fectures) <- colnames(fectures2)
    }
    
    if (method == "SUM") {
        x.s <- matrix(colSums(object, na.rm = TRUE), nrow = 1)
        uni = matrix(rep(1, nrow(object)), ncol = 1)
        area.uni <- uni %*% x.s
        norm.fectures <- (object/area.uni) * 100
    }
    if (method == "NONE") {
        norm.fectures <- object
    }
    
    normtarget <- norm.fectures
    return(normtarget)
}

