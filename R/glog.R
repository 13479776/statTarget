glog <- function(y, lambda) {
    yt <- log(y + sqrt(y^2 + lambda))
    return(yt)
}
