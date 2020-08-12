### Missing value imputation Miss value will be replaced with median value.  degree is the group
### factor mv is th data file A imputed matrix missvalue(mv, degree)
missvalue <- function(mv, degree) {
    for (i in 1:dim(mv)[1]) {
        for (j in 2:dim(mv)[2]) {
            if (mv[i, j] == "0") {
                
                mv[i, j] <- tapply(mv[, j], degree, median, na = TRUE)[degree[i]]
                
            }
        }
    }
    return(mv)
}
