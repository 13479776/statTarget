#' @name pvimPlot
#' @title Gini importance and permutation-based variable 
#' importance measures plots
#' @description Create plots for Gini importance and permutation-based variable 
#' Gini importance measures.
#' @param rForest an object of class randomForest that contains the proximity 
#' component from statTarget rForest function.
#' @param pimpModel an object of permutation-based variable Gini importance 
#' measures (PIMP-algorithm) from statTarget rForest function.
#' @param nvarRF The number of variables in importance plot of randomForest.
#' @param border The color to be used for the border of the bars. 
#' Use border = NA to omit borders. see also barplot.
#' @param space The amount of space (as a fraction of the average bar width) 
#' left before each bar. May be given as a single number or one number per bar.
#' see also barplot
#' @param ...  A generic barplot function from graphics package.
#' @usage pvimPlot(rForest,pimpModel,nvarRF = 6,border= NA,
#' space = 0.3,...)
#' @return The output of the name of selected variable importance.
#' @examples 
#' datpath <- system.file('extdata',package = 'statTarget')
#' statFile <- paste(datpath,'data_example.csv', sep='/')
#' getFile <- read.csv(statFile,header=TRUE)
#' rFtest <- rForest(getFile,ntree = 10,times = 5)
#' pvimPlot(rFtest$randomForest,rFtest$pimpTest)
#' @author Hemi Luan, hemi.luan@gmail.com
#' @export
pvimPlot <- function(rForest, pimpModel, nvarRF = 6, border = NA, space = 0.3, ...) {
    
    gini <- data.frame(rForest$importance[, "MeanDecreaseGini"])
    colnames(gini) <- "GiniImportance"
    gini = cbind(vehicle = row.names(gini), gini)
    if (nvarRF > dim(gini)[1]) {
        stop(paste("The number of variables should not be less than", " nvarRF value (", nvarRF, ")", 
            sep = ""))
    }
    
    sort.gini <- plyr::arrange(gini, plyr::desc(GiniImportance))
    # permutation
    giniPer <- data.frame(pimpModel$pvalue)
    giniPer = cbind(vehicle = row.names(giniPer), giniPer)
    sort.giniP <- data.frame(plyr::arrange(giniPer, p.value))
    sort.giniP[sort.giniP == 0] <- 0.5 * min(sort.giniP[, 2][sort.giniP[, 2] > 0], na.rm = FALSE)
    logtansform <- data.frame(-log10(sort.giniP$p.value))
    sort.giniP.log <- data.frame(cbind(sort.giniP$vehicle, logtansform))
    colnames(sort.giniP.log) <- c("vehicle", "-log10(p.value)")
    # plot
    
    layout(matrix(c(1, 2, 3), 1), c(3, 3, 1))
    # layout(matrix(c(1,2),ncol=2))
    colfunc <- colorRampPalette(c("#009E73", "#D55E00"))(nvarRF)
    
    # barDat_gini <- plyr::arrange(sort.gini[1:nvarRF,2])
    op <- par(mar = c(4.5, 4, 3, 0.5))
    barGini <- graphics::barplot(rev(sort.gini[1:nvarRF, 2]), horiz = TRUE, beside = TRUE, names.arg = rev(sort.gini[1:nvarRF, 
        1]), col = colfunc, space = space, border = border, ylab = "Variables", las = 2, cex.lab = 0.6, 
        cex.main = 0.6, cex.axis = 0.6, cex.names = 0.5, xlab = "Gini Importance", ...)
    text(as.matrix(rev(sort.gini[1:nvarRF, 2])), barGini, labels = round(rev(sort.gini[1:nvarRF, 2]), 
        2), xpd = TRUE, cex = 0.3, col = "grey48", pos = 4)
    par()
    dat <- as.matrix(rev(sort.giniP.log[1:nvarRF, 2]))
    par(op)
    barPvalue <- graphics::barplot(dat, horiz = TRUE, beside = TRUE, space = space, border = border, 
        names.arg = rev(sort.giniP.log[1:nvarRF, 1]), col = colfunc, ylab = "Variables (-log10(p.value))", 
        las = 2, cex.lab = 0.6, cex.main = 0.6, cex.axis = 0.6, cex.names = 0.5, xlab = "Pvalue-based Importance", 
        ...)
    text(dat, barPvalue, labels = round(rev(sort.giniP.log[1:nvarRF, 2]), 2), xpd = TRUE, cex = 0.3, 
        col = "grey48", pos = 4)
    par()
    
    par(mar = c(25, 2, 60, 12), cex = 0.3)
    # par(op1)
    color.bar(colorRampPalette(c("#009E73", "#D55E00"))(6), min = min(sort.giniP.log[1:nvarRF, 2]), 
        max = max(sort.giniP.log[1:nvarRF, 2]))
    return(list(giniImp = sort.gini[1:nvarRF, 1], pvalueImp = sort.giniP.log[1:nvarRF, 1]))
}

pimpRF <- function(X, y, rForest, times = 100, seed = 123, gDist = FALSE) {
    
    ################################################################## randomForest?
    if (!inherits(rForest, "randomForest")) 
        stop("rForest is not of class randomForest")
    mtry = rForest$mtry
    ntree = rForest$ntree
    n = nrow(X)
    p = ncol(X)
    
    set.seed(seed)
    classRF = is.factor(y)
    if (classRF) {
        if (rForest$type == "regression") 
            stop("rForest$type = regression !! y a factor ")
        apply_pb <- function(X, MARGIN, FUN, ...) {
            env <- environment()
            pb_Total <- sum(dim(X)[MARGIN])
            counter <- 0
            pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
            
            wrapper <- function(...) {
                curVal <- get("counter", envir = env)
                assign("counter", curVal + 1, envir = env)
                setTxtProgressBar(get("pb", envir = env), curVal + 1)
                FUN(...)
            }
            res <- apply(X, MARGIN, wrapper, ...)
            close(pb)
            res
        }
        
        y.s = replicate(times, sample(as.integer(y)))
        # Number of permutation
        varip <- apply_pb(y.s, 2, function(y) {
            randomForest::randomForest(X, as.factor(y), mtry = mtry, importance = TRUE, ntree = ntree)[[9]][, 
                "MeanDecreaseGini"]
        })
    } else {
        if (rForest$type == "classification") 
            stop("rForest$type = classification !! y not a factor ")
        y.s = replicate(times, sample(y))
        varip <- apply_pb(y.s, 2, function(y) {
            randomForest::randomForest(X, as.factor(y), mtry = mtry, importance = TRUE, ntree = ntree)[[7]][, 
                "IncNodePurity"]
        })
    }
    # 
    dimNames = dimnames(rForest$importance)[[1]]
    Pimp = list(VarImp = matrix(rForest$importance[, "MeanDecreaseGini"], ncol = 1, dimnames = list(dimNames, 
        "VarImp")), PerVarImp = varip, type = if (classRF) {
        "classification"
    } else {
        "regression"
    }, call = match.call())
    p = nrow(Pimp$PerVarImp)
    if (gDist) {
        mean.PerVarImp = apply(Pimp$PerVarImp, 1, mean)
        sd.PerVarImp = apply(Pimp$PerVarImp, 1, sd)
        j = 1:p
        test.norm = sapply(j, function(j) {
            stats::ks.test(unique(Pimp$PerVarImp[j, ]), "pnorm", mean = mean.PerVarImp[j], sd = sd.PerVarImp[j])$p.value
        })
        
        # the p-value is the probability of observing the original VarImp or one more extreme,
        p.val = sapply(j, function(j) {
            stats::pnorm(Pimp$VarImp[j], mean.PerVarImp[j], sd.PerVarImp[j], lower.tail = FALSE)
        })
        
    } else {
        j = 1:p
        Fn0 = sapply(j, function(j) {
            stats::ecdf(Pimp$PerVarImp[j, ])
        })
        p.val = sapply(j, function(j) {
            1 - Fn0[[j]](Pimp$VarImp[j])
        })
    }
    dimNames = dimnames(Pimp$VarImp)[[1]]
    outTest = list(VarImp = Pimp$VarImp, PerVarImp = Pimp$PerVarImp, gDist = gDist, meanPerVarImp = if (gDist) {
        matrix(mean.PerVarImp, ncol = 1, dimnames = list(dimNames, "mean(PerVarImp)"))
    } else {
        NULL
    }, sdPerVarImp = if (gDist) {
        matrix(sd.PerVarImp, ncol = 1, dimnames = list(dimNames, "sd(PerVarImp)"))
    } else {
        NULL
    }, p.ks.test = if (gDist) {
        matrix(test.norm, ncol = 1, dimnames = list(dimNames, "ks.test"))
    } else {
        NULL
    }, pvalue = matrix(p.val, ncol = 1, dimnames = list(dimNames, "p-value")), type = Pimp$type, call.PIMP = Pimp$call, 
        call = match.call())
    return(outTest)
}

# Function to plot color bar
color.bar <- function(lut, min, max = -min, nticks = 10, ticks = seq(min, max, len = nticks), title = "") {
    scale = (length(lut) - 1)/(max - min)
    plot(c(0, 10), c(min, max), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
        main = title)
    axis(4, round(ticks, 1), las = 1)
    for (i in 1:(length(lut) - 1)) {
        y = (i - 1)/scale + min
        graphics::rect(0, y, 10, y + 1/scale, col = lut[i], border = NA, angle = 0)
    }
}
