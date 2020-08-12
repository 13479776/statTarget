### Description statistic bStatX provide the basic data description for each features.  Data
### descriptions includes mean value, median value, sum, quartile, standard derivatives, etc.  file
### The file with the expression information.  A matrix for data description bStatX(file)
bStatX <- function(file) {
    # xr <- read.csv(file, sep=',', header=TRUE)
    xr = file
    xs = xr[, 3:ncol(xr)]
    x = cbind(xr[, 2], xr[, 1], xs)
    x.nn = x
    sorted = x.nn[order(x.nn[, 1]), ]
    g = c()
    for (i in 1:nrow(sorted)) {
        if (any(g == as.numeric(sorted[i, 1]))) {
            g = g
        } else {
            g = matrix(c(g, as.numeric(sorted[i, 1])), ncol = 1)
        }
    }
    dirout.g = paste(getwd(), "/statTarget/statAnalysis/tmp", sep = "")
    dir.create(dirout.g)
    slink = paste(getwd(), "/statTarget/statAnalysis/DataPretreatment", "/slink.csv", sep = "")
    slink = read.csv(slink, header = TRUE)
    for (i in 1:nrow(g)) {
        vuota <- c()
        fin = matrix(rep(NA, ncol(sorted)), nrow = 1)
        for (j in 1:nrow(sorted)) {
            if (sorted[j, 1] == i) {
                vuota <- as.matrix(sorted[j, ], nrow = 1)
                rownames(vuota) = rownames(sorted)[j]
                fin = rbind(fin, vuota)
            }
        }
        
        nam = paste("r", ExcName(i, slink), sep = ".")
        n = matrix(fin[-1, ], ncol = ncol(sorted))
        n.x = matrix(n[, -1], ncol = ncol(sorted) - 1)
        colnames(n.x) = colnames(x.nn[, 2:ncol(x.nn)])
        name = as.matrix(assign(nam, n.x))
        outputfileg = paste("r.", ExcName(i, slink), ".csv", sep = "")
        write.csv(name, paste(dirout.g, outputfileg, sep = "/"), row.names = FALSE)
    }
    dirout.w = paste(getwd(), "/statTarget/statAnalysis/dataSummary", sep = "")
    dir.create(dirout.w)
    NoF = nrow(g)
    for (i in 1:NoF) {
        ni = paste("r.", ExcName(i, slink), ".csv", sep = "")
        pwdi = paste(getwd(), "/statTarget/statAnalysis/tmp/", ni, sep = "")
        I = read.csv(pwdi, header = TRUE)
        I = I[, -1]
        bS = bStat(I)
        bStat.i = paste("dataSumm_", ExcName(i, slink), ".csv", sep = "")
        assign(bStat.i, bS)
        write.csv(t(bS), paste(dirout.w, bStat.i, sep = "/"))
    }
    tmpfile = paste(getwd(), "/tmp/", sep = "")
    unlink(tmpfile, recursive = TRUE)
}


# This library is free software; you can redistribute it and/or modify it under the terms of the
# GNU Library General Public License as published by the Free Software Foundation; either version
# 2 of the License, or (at your option) any later version.  This library is distributed in the
# hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Library General Public License
# for more details.  You should have received A copy of the GNU Library General Public License
# along with this library; if not, write to the Free Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA


################################################################################ FUNCTION: BASIC STATISTICS: basicStats Returns a basic statistics summary


bStatCor = function(x, ci = 0.95) {
    # A function implemented by Diethelm Wuertz
    
    # Description: Calculates Basic Statistics
    
    # Arguments: x - an object which can be transformed by the function as.matrix() into an object of
    # class matrix.  ci - a numeric value setting the confidence interval.
    
    # Value: a two-column data frame, where the first column takes the value of the statistics, and
    # the second its name, e.g.  'nobs', 'NAs', 'Minimum', 'Maximum', '1. Quartile', '3.  Quartile',
    # 'Mean', 'Median', 'Sum', 'SE Mean', 'LCL Mean', 'UCL Mean', 'Variance', 'Stdev', 'Skewness',
    # 'Kurtosis')
    
    # FUNCTION:
    
    # Univariate/Multivariate:
    y = as.matrix(x)
    
    
    # Handle Column Names:
    if (is.null(colnames(y))) {
        Dim = dim(y)[2]
        if (Dim == 1) {
            colnames(y) = paste(substitute(x), collapse = ".")
        } else if (Dim > 1) {
            colnames(y) = paste(paste(substitute(x), collapse = ""), 1:Dim, sep = "")
        }
    }
    
    # Internal Function - CL Levels:
    cl.vals = function(x, ci) {
        x = x[!is.na(x)]
        n = length(x)
        if (n <= 1) 
            return(c(NA, NA))
        se.mean = sqrt(stats::var(x)/n)
        t.val = qt((1 - ci)/2, n - 1)
        mn = mean(x)
        lcl = mn + se.mean * t.val
        ucl = mn - se.mean * t.val
        c(lcl, ucl)
    }
    
    # Basic Statistics:
    nColumns = dim(y)[2]
    ans = NULL
    for (i in 1:nColumns) {
        X = as.numeric(y[, i])
        # Observations:
        X.length = length(X)
        X = X[!is.na(X)]
        X.na = X.length - length(X)
        # Basic Statistics:
        z = c(X.length, X.na, min(X), max(X), as.numeric(quantile(X, prob = 0.25, na.rm = TRUE)), 
            as.numeric(quantile(X, prob = 0.75, na.rm = TRUE)), mean(X), median(X), sum(X), sqrt(stats::var(X)/length(X)), 
            cl.vals(X, ci)[1], cl.vals(X, ci)[2], stats::var(X), sqrt(stats::var(X)), sqrt(stats::var(X))/mean(X))
        # Row Names:
        znames = c("nobs", "NAs", "Minimum", "Maximum", "1. Quartile", "3. Quartile", "Mean", "Median", 
            "Sum", "SE Mean", "LCL Mean", "UCL Mean", "Variance", "Stdev", "RSD")
        # Output as data.frame
        result = matrix(z, ncol = 1)
        row.names(result) = znames
        ans = cbind(ans, result)
    }
    
    # Column Names:
    colnames(ans) = colnames(y)
    
    # Return Value:
    data.frame(round(ans, digits = 6))
}


################################################################################ 


bStat = function(x, ci = 0.95) {
    # A function implemented by Diethelm Wuertz
    
    # Description: Calculates Basic Statistics
    
    # Arguments: x - an object which can be transformed by the function as.matrix() into an object of
    # class matrix.  ci - a numeric value setting the confidence interval.
    
    # Value: a two-column data frame, where the first column takes the value of the statistics, and
    # the second its name, e.g.  'nobs', 'NAs', 'Minimum', 'Maximum', '1. Quartile', '3.  Quartile',
    # 'Mean', 'Median', 'Sum', 'SE Mean', 'LCL Mean', 'UCL Mean', 'Variance', 'Stdev', 'Skewness',
    # 'Kurtosis')
    
    # FUNCTION:
    
    # Univariate/Multivariate:
    y = as.matrix(x)
    
    
    # Handle Column Names:
    if (is.null(colnames(y))) {
        Dim = dim(y)[2]
        if (Dim == 1) {
            colnames(y) = paste(substitute(x), collapse = ".")
        } else if (Dim > 1) {
            colnames(y) = paste(paste(substitute(x), collapse = ""), 1:Dim, sep = "")
        }
    }
    
    # Internal Function - CL Levels:
    cl.vals = function(x, ci) {
        x = x[!is.na(x)]
        n = length(x)
        if (n <= 1) 
            return(c(NA, NA))
        se.mean = sqrt(stats::var(x)/n)
        t.val = qt((1 - ci)/2, n - 1)
        mn = mean(x)
        lcl = mn + se.mean * t.val
        ucl = mn - se.mean * t.val
        c(lcl, ucl)
    }
    
    # Basic Statistics:
    nColumns = dim(y)[2]
    ans = NULL
    for (i in 1:nColumns) {
        X = y[, i]
        # Observations:
        X.length = length(X)
        X = X[!is.na(X)]
        X.na = X.length - length(X)
        # Basic Statistics:
        z = c(X.length, X.na, min(X), max(X), as.numeric(quantile(X, prob = 0.25, na.rm = TRUE)), 
            as.numeric(quantile(X, prob = 0.75, na.rm = TRUE)), mean(X), median(X), sum(X), sqrt(stats::var(X)/length(X)), 
            cl.vals(X, ci)[1], cl.vals(X, ci)[2], stats::var(X), sqrt(stats::var(X)))
        # Row Names:
        znames = c("nobs", "NAs", "Minimum", "Maximum", "1. Quartile", "3. Quartile", "Mean", "Median", 
            "Sum", "SE Mean", "LCL Mean", "UCL Mean", "Variance", "Stdev")
        # Output as data.frame
        result = matrix(z, ncol = 1)
        row.names(result) = znames
        ans = cbind(ans, result)
    }
    
    # Column Names:
    colnames(ans) = colnames(y)
    
    # Return Value:
    data.frame(round(ans, digits = 6))
}


################################################################################ 









