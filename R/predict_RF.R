#' @name predict_RF
#' @title preidict function for random forest objects in statTarget
#' @description  Prediction of test data using random forest in statTarget.
#' @param object An object created by the function statTarget_rForest.
#' @param newdata A data frame or matrix containing new data. (Note: If not 
#' given, the out-of-bag prediction in object is returned. see randomForest 
#' package.
#' @param type One of response, prob. or votes, indicating the type of output: 
#' predicted values, matrix of class probabilities, or matrix of vote counts. 
#' class is allowed, but automatically converted to 'response', for backward 
#' compatibility.
#' @param ... A generic predict function from randomForest package.
#' @return A class of predicted values is returned. 
#' Object type is classification, for detail see randomForest package.
#' @usage predict_RF(object, newdata, type='response',...)
#' @examples 
#' datpath <- system.file('extdata',package = 'statTarget')
#' statFile <- paste(datpath,'data_example.csv', sep='/')
#' getFile <- read.csv(statFile,header=TRUE)
#' rFtest <- rForest(getFile,ntree = 10,times = 5)
#' predictOutput <- predict_RF(rFtest, getFile[1:19,3:8])
#' @author Hemi Luan, hemi.luan@gmail.com
#' @seealso randomForest
#' @export 
predict_RF <- function(object, newdata, type = "response", ...) {
    # require(randomForest) need object from statTarget_rForest
    if (!inherits(object$randomForest, "randomForest")) 
        stop("Object should be from statTarget_rForest")
    newModel <- object$randomForest
    output <- predict(newModel, newdata, type = "response", ...)
    return(output)
}
