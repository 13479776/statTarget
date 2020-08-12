#' @name rForest
#' @title Random Forest classfication in statTarget
#' @description  rForest provides the Breiman's random 
#' forest algorithm for classification and permutation-based variable 
#'importance measures (PIMP-algorithm).
#' @param file An data frame or 'Stat File' from statTarget software.
#' @param ntree Number of trees to grow. This should not be set to too small a 
#' number, to ensure that every input row gets predicted at least a few times.
#' @param times The number of permutations for 
#' permutation-based variable importance measures. 
#' @param gDist If gDist is TRUE the null importance distributions are 
#' approximated with Gaussian distributions else with empirical cumulative 
#' distributions.
#' @param seed For the same set of random variables and reproducible results.
#' @param ...  A generic function in randomForest package
#' @return Objects  Two objects from statTarget_rForest (1. randomForest,rfModel; 
#' 2. PIMPresult, pimpModel)
#' @return VarImp  The original Gini importance
#' @return PerVarImp  A matrix, where the permuted VarImp measures for the
#' predictor variable.
#' @return p-value  The probability of observing the original VarImp or a larger
#' value, given the fitted null importance distribution.
#' @return 
#' p.ks.test The p-values of the Kolmogorov-Smirnov Tests for 
#' each row PerVarImp. 
#' @usage rForest(file,ntree = 100,times = 100, gDist = TRUE,
#' seed = 123,...)
#' @examples 
#' datpath <- system.file('extdata',package = 'statTarget')
#' statFile <- paste(datpath,'data_example.csv', sep='/')
#' getFile <- read.csv(statFile,header=TRUE)
#' rFtest <- rForest(getFile,ntree = 10,times = 5)
#' @references 
#' Altmann A.,Tolosi L.,Sander O. and Lengauer T. (2010) Permutation importance: 
#' a corrected feature importance measure, Bioinformatics 26 (10), 1340-1347.
#' @references 
#' Ender Celik. (2015) vita: Variable Importance Testing Approaches. R package 
#' version 1.0.0 https://CRAN.R-project.org/package=vita
#' @author Hemi Luan, hemi.luan@gmail.com
#' @export 
rForest <- function(file, ntree = 100, times = 100, gDist = TRUE, seed = 123, ...) {
    x = file
    x.x = x[, 3:ncol(x)]
    rownames(x.x) = x[, 1]
    k = matrix(x[, 2], ncol = 1)
    x.n = cbind(k, x.x)
    sorted = x.n[order(x.n[, 1]), ]
    g = c()
    for (i in 1:nrow(sorted)) {
        if (any(g == sorted[i, 1])) {
            g = g
        } else {
            g = matrix(c(g, sorted[i, 1]), ncol = 1)
        }
    }
    NoF = nrow(g)
    if (NoF >= 2) {
        yF <- as.factor(k[, 1])
        rfModel <- randomForest::randomForest(x.x, yF, ntree = ntree, importance = TRUE, proximity = TRUE, 
            ...)
        pimpModel <- pimpRF(x.x, yF, rfModel, times = times, seed = seed, gDist = gDist)
    } else {
        stop("Single group is not allowed !")
    }
    return(list(randomForest = rfModel, pimpTest = pimpModel))
}

#' @name mdsPlot
#' @title MDSplot in statTarget
#' @description  Multi-dimensional scaling plot of proximity matrix from 
#' randomForest.
#' @param rForest An object of class randomForest that contains the proximity 
#' component from statTarget_rForest function.
#' @param pimpModel An object of permutation-based variable Gini importance 
#' measures (PIMP-algorithm) from statTarget_rForest function.
#' @param Labels Labels is TRUE for visible the sample name in the figure else 
#' with the index for class.
#' @param slink Logical indicating if slinkDat is active for extenal classID.
#' @param slinkDat A data frame for the extenal classID.
#' @param ...  A generic MDSplot function in randomForest package
#' @usage mdsPlot(rForest,pimpModel,Labels = TRUE,slink = FALSE, 
#' slinkDat, ...)
#' @return The output of cmdscale on 1 - rf$proximity is returned invisibly.
#' @seealso MDSplot
#' @examples 
#' datpath <- system.file('extdata',package = 'statTarget')
#' statFile <- paste(datpath,'data_example.csv', sep='/')
#' getFile <- read.csv(statFile,header=TRUE)
#' rFtest <- rForest(getFile,ntree = 10,times = 5)
#' mdsPlot(rFtest$randomForest,rFtest$pimpTest)
#' @author Hemi Luan, hemi.luan@gmail.com
#' @export
mdsPlot <- function(rForest, pimpModel, Labels = TRUE, slink = FALSE, slinkDat, ...) {
    
    mdsData <- stats::cmdscale(1 - rForest$proximity, eig = TRUE, k = 2)
    colnames(mdsData$points) <- paste("Dim", 1:2)
    OOBer <- round(rForest$err.rate[rForest$ntree, "OOB"] * 100, digits = 2)
    if (slink) {
        slink <- slinkDat
    } else {
        xF <- rForest$y
        txF <- factor(xF)
        for (f in 1:length(levels(txF))) {
            levels(txF)[f] <- f
        }
        slink <- data.frame(cbind(xF, as.numeric(txF)))
    }
    if (is.data.frame(slink)) {
        NULL
    } else {
        stop("The slink Data should be set to data frame")
    }
    
    k <- slink[, 2]
    tutticolors = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, "rosybrown4", "green4", "navy", "purple2", "orange", 
        "pink", "chocolate2", "coral3", "khaki3", "thistle", "turquoise3", "palegreen1", "moccasin", 
        "olivedrab3", "azure4", "gold3", "deeppink"), ncol = 1)
    col = c()
    for (s in 1:length(k)) {
        col = c(col, tutticolors[k[s], ])
    }
    mds <- mdsData$points
    lim = c()
    max.pc1 = 1.3 * (max(abs(mds[, 1])))
    max.pc2 = 1.3 * (max(abs(mds[, 2])))
    if (max.pc1 > max.pc2) {
        lim = c(-max.pc1, max.pc1)
    } else {
        lim = c(-max.pc2, max.pc2)
    }
    
    
    Epaxis <- dataEllipse_sT(mds[, 1], mds[, 2], levels = c(0.95), add = FALSE, draw = FALSE, col = "black", 
        lwd = 0.4, plot.points = FALSE, center.cex = 0.2)
    
    if (1.1 * max(Epaxis[, "x"]) < max(lim) & 1.1 * max(Epaxis[, "y"]) < max(lim)) {
        lim = lim
    } else {
        lim = c(1.2 * min(as.numeric(Epaxis)), 1.2 * max(as.numeric(Epaxis)))
    }
    
    # MDS
    mdsPDF <- graphics::plot(mds[, 1], mds[, 2], col = col, pch = 19, xlim = lim, ylim = lim, xlab = "Dim 1", 
        ylab = "Dim 2", sub = paste("OOB estimate of error rate = ", round(OOBer, 2), " %", sep = ""), 
        main = paste("RandomForest MDSPlot", sep = ""), ...)
    axis(1, at = lim * 2, pos = c(0, 0), labels = FALSE, col = "grey", lwd = 0.7)
    axis(2, at = lim * 2, pos = c(0, 0), labels = FALSE, col = "grey", lwd = 0.7)
    dataEllipse_sT(mds[, 1], mds[, 2], levels = c(0.95), add = TRUE, col = "grey48", lwd = 0.4, plot.points = FALSE, 
        center.cex = 0.2)
    
    if (Labels) {
        text(mds[, 1], mds[, 2], mdsPDF, col = col, cex = 0.5, labels = rownames(mds), pos = 1)
    } else {
        legend("topright", legend = levels(factor(slink[, 1])), bty = "n", pch = 19, col = seq_along(levels(factor(slink[, 
            1]))), text.col = seq_along(levels(factor(slink[, 1]))), cex = 0.8)
    }
    return(mds)
}


