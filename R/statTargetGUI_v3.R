#' @name statTargetGUI 
#' @title statTargetGUI for statTarget software
#' @description The statTarget GUI session. The Signal Correction and 
#' Statistical Analysis session are included in statTarget 2.0 software. 
#' See the details at https://stattarget.github.io
#' @author Hemi Luan hemi.luan@gmail.com
#' @references  
#' Dunn WB., et al. Nat Protoc. 2011, 6, pp1060.
#' Luan H., et al. GigaScience 2015, 4, pp16.
#' Luan H., et al. J. Proteome Res., 2015, 14, pp467.
#' @keywords GUI
#' @keywords A GUI of statTarget
#' @keywords statTarget
#' @keywords Metabolomics
#' @keywords Ensemble learning for signal correction
#' @keywords Statistical analysis
#' @return The output of GUI
#' @details  RGTK2 and GTK+ are required for statTargetGUI. We recommend the 64-bit 
#' version of Windows 7, or newer for most users for RGTK2 installation. 
#' For mac OS PC, XQuartz instead of X11 support should be installed for the 
#' Graphical User Interface (GUI). The R 3.3.0 and RGtk2 2.20.31 sailed through the test.
#' @examples 
#' if (interactive()) {statTargetGUI()}
#' @export
statTargetGUI <- function() {
    
    # load_require RGtk2 package
    if (requireNamespace("RGtk2", quietly = TRUE)) {
        print("RGtk2 is loaded correctly")
    } else {
        print("trying to install RGtk2. RGtk2 is preferred on Windows ")
        utils::install.packages("RGtk2")
        if (requireNamespace("RGtk2", quietly = TRUE)) {
            print("RGtk2 installed and loaded")
        } else {
            stop("could not install RGtk2")
        }
    }
    
    
    if (require("gWidgets2RGtk2")) {
        print("gWidgets2RGtk2 is loaded correctly")
    } else {
        print("trying to install gWidgets2RGtk2")
        utils::install.packages("gWidgets2RGtk2")
        if (require("gWidgets2RGtk2")) {
            print("gWidgets2RGtk2 installed and loaded")
        } else {
            stop("could not install gWidgets2RGtk2")
        }
    }
    
    if (require("gWidgets2")) {
        print("gWidgets2 is loaded correctly")
    } else {
        print("trying to install gWidgets2")
        utils::install.packages("gWidgets2")
        if (require("gWidgets2")) {
            print("gWidgets2 installed and loaded")
        } else {
            stop("could not install gWidgets2")
        }
    }
    
    # GUI
    
    color <- grDevices::colors()
    linetype <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F8", "431313", "22848222")
    linetype <- rep(linetype, 4)
    widgets <- list()
    win = gWidgets2::gwindow("Welcome to statTarget", expand = FALSE, fill = FALSE, height = 500, width=250)
    
    gp = gWidgets2::ggroup(horizontal = FALSE, container = win, expand = FALSE, fill = FALSE)
    gi = gWidgets2::gimage("shinv.png", system.file("extdata", package = "statTarget"), cont = gp)
    
    tmp <- gWidgets2::gframe("Data-inputs", container = gp, expand = FALSE, fill = TRUE)
    lyout <- gWidgets2::glayout(container = tmp)
    lyout[1, 1] <- gWidgets2::gbutton("Meta File...", cont = lyout, handler = function(h, ...) {
        std <- gWidgets2::gfile(text = "Select Meta-information File...", filter = list(metafiles = list(patterns = c("*.csv")), 
            `All files` = list(patterns = c("*"))))
        if (std != "") {
            if (length(grep("\\", std, fixed = TRUE)) > 0) {
                mystr <- strsplit(std, split = "\\", fixed = TRUE)[[1]]
                mystr.lth <- mystr[length(mystr)]
                mydir <- substr(std, 1, stop = (nchar(std) - nchar(mystr.lth) - 1))
                setwd(mydir)
            } else if (length(grep("/", std, fixed = TRUE)) > 0) {
                mystr <- strsplit(std, split = "/", fixed = TRUE)[[1]]
                mystr.lth <- mystr[length(mystr)]
                mydir <- substr(std, 1, stop = (nchar(std) - nchar(mystr.lth) - 1))
                setwd(mydir)
            }
            samP <- std
            gWidgets2::svalue(widgets$pheno) <- samP
        }
    })
    lyout[1, 2] <- (widgets$pheno <- gWidgets2::gedit(text = "", initial.msg = "Click the left button", 
        cont = lyout))
    gWidgets2::tooltip(widgets$pheno) <- "Meta-information for Signal Correction, i.e. sample name, class, batch and order"
    
    lyout[2, 1] <- gWidgets2::gbutton("Profile File...", cont = lyout, handler = function(h, ...) {
        std <- gWidgets2::gfile("Select Profile File...", filter = list(`Profile files` = list(patterns = c("*.csv")), 
            `All files` = list(patterns = c("*"))))
        if (std != "") {
            if (length(grep("\\", std, fixed = TRUE)) > 0) {
                mystr <- strsplit(std, split = "\\", fixed = TRUE)[[1]]
                mystr.lth <- mystr[length(mystr)]
                mydir <- substr(std, 1, stop = (nchar(std) - nchar(mystr.lth) - 1))
                setwd(mydir)
            } else if (length(grep("/", std, fixed = TRUE)) > 0) {
                mystr <- strsplit(std, split = "/", fixed = TRUE)[[1]]
                mystr.lth <- mystr[length(mystr)]
                mydir <- substr(std, 1, stop = (nchar(std) - nchar(mystr.lth) - 1))
                setwd(mydir)
            }
            samF <- std
            gWidgets2::svalue(widgets$profile) <- samF
        }
    })
    lyout[2, 2] <- (widgets$profile <- gWidgets2::gedit(text = "", initial.msg = "Click the left button", 
        cont = lyout))
    gWidgets2::tooltip(widgets$profile) <- "Expression data for Signal Correction"
    
    
    
    ## prepare for data file ready
    nb = gWidgets2::gnotebook(cont = gp, fill = TRUE, expand = TRUE)
    shiftco_win = gWidgets2::ggroup(horizontal = FALSE, cont = nb, label = "Signal Correction")
    stat_win = gWidgets2::ggroup(horizontal = FALSE, cont = nb, label = "Statistical Analysis")
    mix = gWidgets2::ggroup(horizontal = FALSE, cont = nb, label = "Trans-->X")
    sb <- gWidgets2::gstatusbar("Contact Us: luanhm@sustc.edu.cn\nSouthern University of Science and Technology", container = win)
    gWidgets2::font(sb) <- list(size = 9, color = "blue")
    svalue(nb) <- 1
    
    # tranx
    gWidgets2::glayout(container = mix)
    xlyout <- gexpandgroup("transX Function", cont = mix)
    tlyout<-gWidgets2::glayout(container=xlyout,spacing = 10)

    tlyout[1,1]<-gWidgets2::gbutton("Data File", 
                                   cont=tlyout, handler = function(h,...) {
                                     stdt<-gWidgets2::gfile(text="Select Data File...",
                                                           filter=list("Data files" = list(patterns = c("*.csv","*.tsv")),
                                                                       "All files"=list(patterns=c("*"))))
                                     if(stdt != "")
                                     {
                                       if(length(grep("\\",stdt,fixed=TRUE))>0){
                                         mystr<-strsplit(stdt,split="\\",fixed=TRUE)[[1]]
                                         mystr.lth<-mystr[length(mystr)]
                                         mydir<-substr(stdt,1,stop=(nchar(stdt)-nchar(mystr.lth)-1))
                                         setwd(mydir)
                                       }else if(length(grep("/",stdt,fixed=TRUE))>0){
                                         mystr<-strsplit(stdt,split="/",fixed=TRUE)[[1]]
                                         mystr.lth<-mystr[length(mystr)]
                                         mydir<-substr(stdt,1,stop=(nchar(stdt)-nchar(mystr.lth)-1))
                                         setwd(mydir)
                                       }
                                       inputdat <- stdt
                                       gWidgets2::svalue(widgets$inputdat) <- inputdat
                                     }
                                   })
    tlyout[1,2]<-(widgets$inputdat<-gWidgets2::gedit(text="",initial.msg = "Click the left button",cont=tlyout))
    gWidgets2::tooltip(widgets$inputdat) <- "The output-file from XCMS, MZmine2,Skyline or SIEVE. "
    tlyout[2,1] <- "DataSource"
    tlyout[2,2] <- widgets$dsour <- gWidgets2::gradio(c("XCMS","MZmine2","SKYLINE","SIEVE"), selected = 1,cont = tlyout, horizontal = FALSE)
    tlyout[2,3] <- button.group <- gWidgets2::ggroup(container = tlyout)
    
    gWidgets2::addSpring(button.group)
    gWidgets2::gbutton("Run", handler=function(h,...){
    transdata <- gWidgets2::svalue(widgets$inputdat)
    typeT <- gWidgets2::svalue(widgets$dsour)
    transX(data=transdata, type= typeT)
    galert(paste(" Task Finished ! See data at\n" 
                 ,getwd()), title = typeT, delay = 5)
    },container=button.group
    )
    gWidgets2::visible(xlyout) <- TRUE
    
    # shiftco_win
    lyout <- gWidgets2::glayout(container = shiftco_win)
    
    # Data preprocessing
    gWidgets2::glayout(container = shiftco_win)
    slyout <- ggroup(cont = shiftco_win, horizontal = FALSE, spacing = 15)
    
    dslyout <- gexpandgroup("Data preprocessing", cont = slyout)
    dslyout2 <- gformlayout(cont = dslyout)
    
    
    widgets$Frule1 <- gWidgets2::gedit("0.8", initial = "0.8", label = "NA.Filter", cont = dslyout2)
    gWidgets2::tooltip(widgets$Frule1) <- "Removing missing values using 80 percent rule; valueRange,'0 ~ 1' "
    
    widgets$imputeM1 <- gWidgets2::gcombobox(c("KNN", "min", "minHalf", "median"), label = "Imputation", 
        cont = dslyout2)
    gWidgets2::tooltip(widgets$imputeM1) <- "Missing value or zero value imputation"
    
    widgets$Plot <- gWidgets2::gcombobox(c("TRUE", "FALSE"), label = "plotQC", 
                                             cont = dslyout2)
    gWidgets2::tooltip(widgets$Plot) <- "Plot the intesity vs injection order of all sampels"
    
    gWidgets2::visible(dslyout) <- TRUE
    
    # Singal correction method
    
    cslyout <- gexpandgroup("Singal correction methods", cont = slyout)
    cslyout2 <- gformlayout(cont = cslyout)
    
    
    widgets$MLmethod <- gWidgets2::gradio(c("QCRFSC", "QCRLSC","Combat"), 
                                          selected = 1, cont = cslyout2, horizontal = FALSE)
    
    gWidgets2::visible(cslyout) <- TRUE
    
    # Parameters for singal correction method
    qfslyout <- gexpandgroup("Parameters", cont = slyout)
    qfslyout2 <- glayout(cont = qfslyout, homogeneous = TRUE)
    
    
    qfslyout2[1, 1, expand = TRUE, anchor = c(1, 0)] <- "Ntree (QCRFSC only)"
    qfslyout2[1, 2] <- (widgets$ntree1 <- gWidgets2::gslider(from = 100, to = 1000, by = 10, value = 500, 
        cont = qfslyout2))
    
    qfslyout2[2, 1, expand = TRUE, anchor = c(1, 0)] <- "      QCspan (QCRLSC only)"
    qfslyout2[2, 2] <- (widgets$QCspan <- gWidgets2::gslider(from = 0, to = 0.75, by = 0.01, value = 0, 
        cont = qfslyout2))
    qfslyout2[3, 1, expand = TRUE, anchor = c(1, 0)] <- "      CV% Cutoff"
    qfslyout2[3, 2] <- (widgets$coCV1 <- gWidgets2::gslider(from = 0, to = 100, by = 5, value = 30, 
                                                             cont = qfslyout2))
    
    
    gWidgets2::visible(qfslyout) <- FALSE
    
    glabel("Click above to show", cont = slyout)
    
    button.group <- gWidgets2::ggroup(container = shiftco_win)
    ## Push buttons to right
    gWidgets2::addSpring(button.group)
    gsc <- gWidgets2::gbutton("Run", handler = function(h, ...) {
        # close.cur.dev()
        samPeno = gWidgets2::svalue(widgets$pheno)
        samFile = gWidgets2::svalue(widgets$profile)
        Frule1 = gWidgets2::svalue(widgets$Frule1)
        Frule1 = as.numeric(Frule1)
        QCspan = gWidgets2::svalue(widgets$QCspan)
        QCspan = as.numeric(QCspan)
        ntree1 = gWidgets2::svalue(widgets$ntree1)
        ntree1 = as.numeric(ntree1)
        MLmethod = gWidgets2::svalue(widgets$MLmethod)
        imputeM1 = gWidgets2::svalue(widgets$imputeM1)
        coCV1 = gWidgets2::svalue(widgets$coCV1)
        coCV1 = as.numeric(coCV1)
        
        Plot = gWidgets2::svalue(widgets$Plot)
        
        # QC span check
        if (MLmethod == "QCRLSC" & as.numeric(svalue(widgets$QCspan)) > 1e-04 & as.numeric(svalue(widgets$QCspan)) < 
            0.3499) {
            
            QCok <- (gWidgets2::gconfirm("ERROR! For avoiding overfitting of the observed data, the default value of QCspan should be set as '0'.  
                     STOP program now ?", 
                icon = c("error"), title = "QCspan Alert"))
            if (QCok) {
                stop("Reset your QCspan please! ")
            } else {
                NULL
            }
        } else {
            NULL
        }
        
        logw1 <- gWidgets2::gwindow("Signal Correction Start...")
        size(logw1) <- c(600, 500)
        logg1 <- gWidgets2::gvbox(container = logw1)
        logg1$set_borderwidth(5)
        
        if(MLmethod %in% c("QCRLSC","QCRFSC")) {
        utils::capture.output(shiftCor(samPeno, samFile, Frule = Frule1, MLmethod = MLmethod, ntree = ntree1, 
            QCspan = QCspan, degree = 2, imputeM = imputeM1, coCV = coCV1,plot = Plot), file = "shiftCor.log", split = TRUE, 
            append = FALSE)
        }
        if(MLmethod %in% c("Combat")) {
          utils::capture.output(shiftCor_dQC(samPeno,samFile, Frule = 0.8, MLmethod = "Combat"),
                                file = "shiftCor.log", split = TRUE, append = FALSE)
        }
        
        logtmp1 <- try(readLines(paste(getwd(), "shiftCor.log", sep = "/")), silent=TRUE)
        if ("try-error" %in% attr(logtmp1,"class")){
          logtmp1 <- "

                       Cannot load the log files!\n 
                       Check your results, please!"
          ltxt1 <- gWidgets2::gtext(paste(logtmpP1 = logtmp1, collapse = "\n"), wrap = FALSE,
                                    cont = logg1, expand = TRUE, fill = TRUE)
          
        } else {
        ltxt1 <- gWidgets2::gtext(paste(logtmpP1 = logtmp1, collapse = "\n"), wrap = FALSE,
                                  cont = logg1, expand = TRUE, fill = TRUE)}
        gWidgets2::gseparator(cont = logg1)
        bglog1 <- gWidgets2::ggroup(cont = logg1)
        gWidgets2::addSpring(bglog1)
        gWidgets2::gbutton("All Done !", cont = bglog1, handler = function(h, ...) {
            dispose(logw1)
        })
        
        # 
        about <- "The shifCor function provides the QC-based signal correction 
        for your omic data. The path of the shifCor output files 
        was showed in this log file. 
        
        See the manual, vignettes and references at 
        https: //stattarget.github.io"
        
        gWidgets2::gbutton("about", cont = bglog1, handler = function(...) {
            w1 <- gwindow("about", parent = logw1, expand = TRUE, fill = FALSE)
            g <- gvbox(cont = w1, space = 10L)
            g$set_borderwidth(5L)
            glabel(about, cont = g, editable = FALSE)
            gseparator(cont = g)
            bg <- ggroup(cont = g)
            addSpring(bg)
            gbutton("dismiss", cont = bg, handler = function(...) dispose(w1))
        })
    }, container = button.group)
    
    
    
    # Data preprocessing
    gWidgets2::glayout(container = stat_win)
    lyoutStat <- gWidgets2::ggroup(cont = stat_win, horizontal = FALSE, spacing = 15)
    sflyout <- gexpandgroup("Data-inputs", cont = lyoutStat)
    sflyout2 <-gWidgets2::glayout(container=sflyout,spacing = 10)
    
    
    
    sflyout2[1, 1] <- gWidgets2::gbutton("Stat File...", cont = sflyout, handler = function(h, ...) {
      std <- gWidgets2::gfile("Select Stat File...", filter = list(`Stat files` = list(patterns = c("*.csv")), 
                                                                   `All files` = list(patterns = c("*"))))
      if (std != "") {
        if (length(grep("\\", std, fixed = TRUE)) > 0) {
          mystr <- strsplit(std, split = "\\", fixed = TRUE)[[1]]
          mystr.lth <- mystr[length(mystr)]
          mydir <- substr(std, 1, stop = (nchar(std) - nchar(mystr.lth) - 1))
          setwd(mydir)
        } else if (length(grep("/", std, fixed = TRUE)) > 0) {
          mystr <- strsplit(std, split = "/", fixed = TRUE)[[1]]
          mystr.lth <- mystr[length(mystr)]
          mydir <- substr(std, 1, stop = (nchar(std) - nchar(mystr.lth) - 1))
          setwd(mydir)
        }
        stat <- std
        gWidgets2::svalue(widgets$stat) <- stat
      }
    })
    sflyout2[1, 2] <- (widgets$stat <- gWidgets2::gedit(text = "", initial.msg = "Click the left button", 
                                                     cont = sflyout))
    gWidgets2::tooltip(widgets$stat) <- "Expression data for Statistical Analysis"
    
    gWidgets2::visible(sflyout) <- TRUE
    
    #lyoutStat <- gWidgets2::ggroup(cont = stat_win, horizontal = FALSE, spacing = 15)
    
    glyoutStat <- gWidgets2::gexpandgroup("Data preprocessing", cont = lyoutStat)
    glyoutStat2 <- gWidgets2::gformlayout(cont = glyoutStat)
    
    widgets$Frule2 <- gWidgets2::gedit("0.8", initial = "0.8", label = "   NA.Filter", cont = glyoutStat2)
    gWidgets2::tooltip(widgets$Frule2) <- "Removing missing values using 80 percent rule; valueRange,'0 ~ 1' "
    
    widgets$imputeM2 <- gWidgets2::gcombobox(c("KNN", "min", "minHalf", "median"), label = "   Imputation", 
        cont = glyoutStat2)
    gWidgets2::tooltip(widgets$imputeM2) <- "Missing value imputation"
    
    widgets$normM <- gWidgets2::gcombobox(c("NONE", "SUM", "PQN"), label = "   Normalization", cont = glyoutStat2)
    gWidgets2::tooltip(widgets$normM) <- "median quotient normalization, 'PQN'; integral normalization , 'SUM', and 'NONE' "
    
    widgets$Glog <- gWidgets2::gcombobox(c("TRUE", "FALSE"), label = "   Glog", cont = glyoutStat2)
    gWidgets2::tooltip(widgets$Glog) <- "Data variance stabilising transformations"
    
    
    gWidgets2::visible(glyoutStat) <- FALSE
    
    
    # PCA and PLS analysis
    plyoutStat <- gWidgets2::gexpandgroup("PCA and PLS analysis", cont = lyoutStat, expand = TRUE)
    plyoutStat2 <- gWidgets2::gformlayout(container = plyoutStat, expand = TRUE)
    
    widgets$scalingMethod <- gWidgets2::gcombobox(c("Center", "Pareto", "Auto", "Vast", "Range", "None"), 
        label = "   Scaling", cont = plyoutStat2)
    gWidgets2::tooltip(widgets$scalingMethod) <- "Data scaling for PCA or PLS(-DA) analysis"
    
    
    widgets$pcaX <- gWidgets2::gedit("1", width = 8, label = "   PCs in Xaxis", cont = plyoutStat2)
    gWidgets2::tooltip(widgets$pcaX) <- "   The X-axis (horizontal) component"
    
    widgets$pcaY <- gWidgets2::gedit("2", width = 8, label = "   PCs in Xaxis", cont = plyoutStat2)
    gWidgets2::tooltip(widgets$pcaY) <- "   The Y-axis (Vertical) component"
    
    gWidgets2::visible(plyoutStat) <- FALSE
    
    # random forest analysis
    rlyoutStat <- gexpandgroup("Random forest analysis", cont = lyoutStat)
    rlyoutStat2 <- gWidgets2::gformlayout(container = rlyoutStat)
    
    widgets$nvarRF <- gWidgets2::gedit("5", width = 8, label = "   nvarRF", cont = rlyoutStat2)
    gWidgets2::tooltip(widgets$nvarRF) <- "Visualizing variables importance in randomforest model"
    
    widgets$ntree2 <- gWidgets2::gedit("500", width = 8, label = "   ntree", cont = rlyoutStat2)
    gWidgets2::tooltip(widgets$ntree2) <- "The number of trees to grow in random forest model"
    
    gWidgets2::visible(rlyoutStat) <- FALSE
    
    # Permutation times
    tlyoutStat <- gWidgets2::gexpandgroup("Permutation", cont = lyoutStat)
    tlyoutStat2 <- gWidgets2::gformlayout(container = tlyoutStat)
    
    widgets$Permutation <- gWidgets2::gedit("20", width = 8, label = "   Permutation times", cont = tlyoutStat2)
    gWidgets2::tooltip(widgets$Permutation) <- "The number of permutation for binary PLS(-DA) and randomForest model"
    
    gWidgets2::visible(tlyoutStat) <- FALSE
    
    # Univariate analysis
    ulyoutStat <- gWidgets2::gexpandgroup("Univariate analysis", cont = lyoutStat)
    ulyoutStat2 <- gWidgets2::glayout(container = ulyoutStat, expand = TRUE, fill = TRUE)
    
    ulyoutStat2[1, 1] <- "   Volcano FC >"
    ulyoutStat2[1, 2] <- (widgets$mfc <- gWidgets2::gedit("2", width = 5, label = "Volcano FC >", 
        cont = ulyoutStat2))
    gWidgets2::tooltip(widgets$mfc) <- "Fold changes threshold for Volcano plot"
    
    ulyoutStat2[2, 1] <- "   Volcano FC <"
    ulyoutStat2[2, 2] <- (widgets$lfc <- gWidgets2::gedit("0.5", width = 5, label = "Volcano FC <", 
        cont = ulyoutStat2))
    gWidgets2::tooltip(widgets$lfc) <- "Fold changes threshold for Volcano plot"
    
    ulyoutStat2[3, 1] <- "   Volcano Pvalue <"
    ulyoutStat2[3, 2] <- (widgets$pvalue <- gWidgets2::gedit("0.05", width = 5, label = "Volcano Pvalue <", 
        cont = ulyoutStat2))
    gWidgets2::tooltip(widgets$pvalue) <- "Significance level 'p-value' for Volcano plot"
    
    ulyoutStat2[1, 3] <- "   Multiple Testing"
    ulyoutStat2[1, 4, expand = TRUE] <- (widgets$FDR <- gWidgets2::gcombobox(c("TRUE", "FALSE"), cont = ulyoutStat2))
    gWidgets2::tooltip(widgets$FDR) <- "Controlling the false discovery rate with Benjamini Hochberg method"
    
    ulyoutStat2[2, 3] <- "   Labels"
    ulyoutStat2[2, 4, expand = TRUE] <- (widgets$Labels <- gWidgets2::gcombobox(c("TRUE", "FALSE"), 
        cont = ulyoutStat2))
    gWidgets2::tooltip(widgets$Labels) <- "Labelling sampleNames for scatterPlots of multiple statistical analysis"
    
    gWidgets2::visible(ulyoutStat) <- FALSE
    
    gWidgets2::glabel("Click above to show", cont = lyoutStat)
    
    button.group <- gWidgets2::ggroup(container = stat_win)
    ## Push buttons to right
    gWidgets2::addSpring(button.group)
    gWidgets2::gbutton("Run", handler = function(h, ...) {
        
        logw2 <- gWidgets2::gwindow("Statistical Analysis Start...")
        size(logw2) <- c(600, 500)
        logg2 <- gWidgets2::gvbox(container = logw2)
        logg2$set_borderwidth(5)
        
        
        # run the code
        file = gWidgets2::svalue(widgets$stat)
        Frule2 = gWidgets2::svalue(widgets$Frule2)
        Frule2 = as.numeric(Frule2)
        imputeM2 = gWidgets2::svalue(widgets$imputeM2)
        glog = gWidgets2::svalue(widgets$Glog)
        normM = gWidgets2::svalue(widgets$norm)
        nvarRF = gWidgets2::svalue(widgets$nvarRF)
        nvarRF = as.numeric(nvarRF)
        ntree2 = gWidgets2::svalue(widgets$ntree2)
        ntree2 = as.numeric(ntree2)
        scaling = gWidgets2::svalue(widgets$scalingMethod)
        silt = gWidgets2::svalue(widgets$Permutation)
        silt = as.numeric(silt)
        pcax = gWidgets2::svalue(widgets$pcaX)
        pcax = as.numeric(pcax)
        pcay = gWidgets2::svalue(widgets$pcaY)
        pcay = as.numeric(pcay)
        Labels = gWidgets2::svalue(widgets$Labels)
        FDR = gWidgets2::svalue(widgets$FDR)
        upper.lim = gWidgets2::svalue(widgets$mfc)
        upper.lim = as.numeric(upper.lim)
        lower.lim = gWidgets2::svalue(widgets$lfc)
        lower.lim = as.numeric(lower.lim)
        sig.lim = gWidgets2::svalue(widgets$pvalue)
        sig.lim = as.numeric(sig.lim)
        
        utils::capture.output(statAnalysis(file = file, Frule = Frule2, normM = normM, imputeM = imputeM2, 
            glog = glog, FDR = FDR, ntree = ntree2, nvarRF = nvarRF, scaling = scaling, silt = silt, 
            pcax = pcax, pcay = pcay, Labels = Labels, upper.lim = upper.lim, lower.lim = lower.lim, 
            sig.lim = sig.lim), file = "statAnalysis.log", split = TRUE, append = FALSE)
        
        
        logtmp2 <- try(readLines(paste(getwd(), "statAnalysis.log", sep = "/")), silent=TRUE)
        if ("try-error" %in% attr(logtmp2,"class")){
          logtmp2 <- "
                      
                      Cannot load the log files!\n 
                      Check your results, please!"
          ltxt <- gtext(paste(logtmpP = logtmp2, collapse = "\n"), wrap = FALSE, 
                        cont = logg2, expand = TRUE, fill = TRUE)
        } else {
          ltxt <- gtext(paste(logtmpP = logtmp2, collapse = "\n"), wrap = FALSE, 
                        cont = logg2, expand = TRUE, fill = TRUE)
        }
        
        gseparator(cont = logg2)
        
        bglog <- ggroup(cont = logg2)
        addSpring(bglog)
        gbutton("All Done !", cont = bglog, handler = function(h, ...) {
            dispose(logw2)
        })
        
        # 
        about <- "The statAnalysis function provides the statistical analysis 
for your omic data. The path of the statAnalsis output files 
was showed in this log file. 

See the manual, vignettes and references at 
https: //stattarget.github.io"
        
        gbutton("about", cont = bglog, handler = function(...) {
            w1 <- gwindow("about", parent = logw2, expand = TRUE, fill = FALSE)
            g <- gvbox(cont = w1, space = 10L)
            g$set_borderwidth(5L)
            glabel(about, cont = g, editable = FALSE)
            gseparator(cont = g)
            bg <- ggroup(cont = g)
            addSpring(bg)
            gbutton("dismiss", cont = bg, handler = function(...) dispose(w1))
        })
        
        
    }, container = button.group)
}
