sT_univariate <- function(file, FDR = FDR, plot.volcano, upper.lim, lower.lim, sig.lim, save.boxplot) {
    dirout.uni = paste(getwd(), "/Univariate/", sep = "")
    dir.create(dirout.uni)
    comp = read.csv(file, sep = ",", header = TRUE)
    comp.x = comp[, 3:ncol(comp)]
    comp.x = cbind(comp[, 2], comp[, 1], comp.x)
    pwdfile = paste(getwd(), "/Univariate/DataTable.csv", sep = "")
    write.csv(comp.x, pwdfile, row.names = FALSE)
    
    # cal min no. of levels
    checkNum <- min(summary(as.factor(comp.x[,1])))
    
    cat("\n", "P-value Calculating...")
    
    if (FDR) {
        cat("\n", "P-value was adjusted using Benjamini-Hochberg Method")
    }
    
    shapiro(pwdfile)
    welch(pwdfile)
    wmw(pwdfile)
    cat("\n")
    cat("\n", "Odd.Ratio Calculating...", "\n")
    oddRatio(pwdfile)
    
    if(checkNum <= 5) {
      cat("\n", "*ROC analysis skipped", "\n")
    } else {
    cat("\n", "ROC Calculating...", "\n")
    aucROC(pwdfile)
    }
    # cat('\n','RandomForest Calculating...','\n')
    
    # RandomF(file,nvarRF)
    
    
    if (FDR) {
        pvalues(pwdfile, fdr = TRUE)
    } else {
        pvalues(pwdfile, fdr = FALSE)
    }
    col_pvalues(pwdfile)
    if (plot.volcano) {
        cat("\n", "Volcano Plot Output...", "\n")
        volcano(pwdfile, upper.lim = upper.lim, lower.lim = lower.lim, sig.lim = sig.lim)
    }
    if (save.boxplot) {
        cat("\n", "Box Plot Output...", "\n")
        boxPlot(pwdfile)
    }
}
