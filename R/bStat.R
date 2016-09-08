bStat = 
  function(x, ci = 0.95) 
  {   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Calculates Basic Statistics
    
    # Arguments:
    #   x - an object which can be transformed by the function
    #       as.matrix() into an object of class matrix. 
    #   ci - a numeric value setting the confidence interval.
    
    # Value:
    #   a two-column data frame, where the first column takes the 
    #   value of the statistics, and the second its name, e.g.
    #   "nobs", "NAs",  "Minimum", "Maximum", "1. Quartile",  
    #   "3. Quartile",  "Mean", "Median", "Sum",  "SE Mean", 
    #   "LCL Mean", "UCL Mean", "Variance", "Stdev", "Skewness", 
    #   "Kurtosis")
    
    # FUNCTION:
    
    # Univariate/Multivariate:
    y = as.matrix(x)
    
    
    # Handle Column Names:
    if (is.null(colnames(y))) {
      Dim = dim(y)[2]
      if (Dim == 1) {
        colnames(y) = paste(substitute(x), collapse = ".")
      } else if (Dim > 1) {
        colnames(y) = 
          paste(paste(substitute(x), collapse = ""), 1:Dim, sep = "")
      }
    }
    
    # Internal Function - CL Levels:    
    cl.vals = function(x, ci) {
      x = x[!is.na(x)]
      n = length(x)
      if(n <= 1) return(c(NA, NA))
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
      z = c(
        X.length, X.na, min(X), max(X),
        as.numeric(quantile(X, prob = 0.25, na.rm = TRUE)), 
        as.numeric(quantile(X, prob = 0.75, na.rm = TRUE)), 
        mean(X), median(X), sum(X), sqrt(stats::var(X)/length(X)), 
        cl.vals(X, ci)[1], cl.vals(X, ci)[2], stats::var(X), 
        sqrt(stats::var(X)))    
      # Row Names:
      znames = c(
        "nobs", "NAs",  "Minimum", "Maximum", 
        "1. Quartile",  "3. Quartile",  "Mean", "Median", 
        "Sum",  "SE Mean", "LCL Mean", "UCL Mean", 
        "Variance", "Stdev")   
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

