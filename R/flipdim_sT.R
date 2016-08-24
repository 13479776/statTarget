flipdim_sT <- function(a, dim=1) {
  if (!is.matrix(a))
    stop("Argument 'a' must ba a matrix.")
  if (!(dim %in% c(1,2)))
    stop("Argument 'dim' must be 1 or 2 (for rows or columns).")
  
  if (dim == 1) {
    a <- a[nrow(a):1, ]
  } else {
    a <- a[ ,ncol(a):1]
  }
  return(a)
}
