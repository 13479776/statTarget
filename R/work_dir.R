#' @importFrom utils read.csv
#' @importFrom utils read.delim
#' @importFrom utils write.csv
#' @importFrom utils write.table
#' @importFrom stats binomial
#' @importFrom utils capture.output
#' @importFrom stats cor
#' @importFrom stats var
## @importFrom stats cov
#' @importFrom stats density
#' @importFrom stats ecdf
#' @importFrom stats glm
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats pf
#' @importFrom stats prcomp
#' @importFrom stats qt
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats shapiro.test
#' @importFrom stats t.test
## @importFrom stats var
#' @importFrom stats wilcox.test
#' @importFrom stats loess
#' @importFrom stats predict
#' @importFrom stats update
#' @importFrom stats resid
#' @importFrom graphics axis
#' @importFrom graphics barplot
#' @importFrom graphics boxplot
#' @importFrom graphics grid
#' @importFrom graphics hist
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics segments
#' @importFrom graphics text
#' @importFrom graphics strwidth
#' @importFrom graphics strheight
#' @importFrom graphics box
#' @importFrom graphics identify
#' @importFrom graphics lines
#' @importFrom graphics title
#' @importFrom graphics mtext
#' @importFrom utils install.packages
#' @importFrom randomForest randomForest
#' @importFrom randomForest MDSplot
#' @importFrom randomForest importance
#' @importFrom grDevices rainbow
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rgb
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.copy2pdf
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.new
#' @importFrom grDevices colors
#' @importFrom grDevices palette
#' @importFrom grDevices col2rgb
# @importFrom gWidgets2 svalue
#' @importFrom plyr rlply
#' @importFrom plyr daply
#' @importFrom plyr .
# @importFrom gWidgets2 gwindow @importFrom gWidgets2 ggroup @importFrom gWidgets2 gimage
# @importFrom gWidgets2 gframe @importFrom gWidgets2 glayout @importFrom gWidgets2 gbutton
# @importFrom gWidgets2 svalue @importFrom gWidgets2 gfile @importFrom gWidgets2 gedit @importFrom
# gWidgets2 font @importFrom gWidgets2 gstatusbar
#' @import pls
#' @import pdist 
#' @importFrom ROC rocdemo.sca
#' @importFrom ROC AUC
#' @importFrom rrcov plot
#' @import impute
## @import stats @import gWidgets2RGtk2
#' @importFrom graphics abline
#' @importFrom stats lm
#' @importFrom stats complete.cases
#' @importFrom stats cov.wt
#' @importFrom stats qf
#' @importFrom graphics layout
#' @importFrom stats aggregate
#' @importFrom utils setTxtProgressBar 
#' @importFrom utils txtProgressBar
work_dir <- function(dir.name) {
    WorkinDir = paste(getwd(), "/", dir.name, "/", sep = "")
    dir.create(WorkinDir)
    file = list.files()
    file.copy(file, WorkinDir)
    setwd(WorkinDir)
    unlink("statTarget.R")
}

