ST_rForest <- function(file, ntree = 500, times = 100, gDist = TRUE, seed = 123, Labels = TRUE, nvarRF = 6) {
    
    dirout.w = paste(getwd(), "/randomForest", sep = "")
    # dirout.RF = paste(getwd(), '/RF_Plots/', sep='')
    dir.create(dirout.w)
    # dir.create(dirout.RF)
    cat("\n")
    statRF <- rForest(file, ntree = ntree, times = times, gDist = gDist, seed = seed)
    
    slinkPath = paste(getwd(), "/DataPretreatment", "/slink.csv", sep = "")
    getSlink = read.csv(slinkPath, header = TRUE)
    
    cat("\n", "Random Forest Model")
    cat("\n", "Type: ", statRF$randomForest$type)
    
    # cat('\n','Response Vector:','\n') pinfo <- unique(getSlink) rownames(pinfo) <- pinfo[,2];
    # pinfo[,2] <- NULL print(t(pinfo))
    cat("\n", "ntree:", statRF$randomForest$ntree)
    cat("\n", "mtry: ", statRF$randomForest$mtry)
    OOB <- round(statRF$randomForest$err.rate[statRF$randomForest$ntree, "OOB"] * 100, digits = 2)
    cat("\n", "OOB estimate of error rate:", OOB, "%", "\n")
    # print(statRF$randomForest$confusion)
    machGini <- data.frame(statRF$randomForest$importance[, "MeanDecreaseGini"])
    colnames(machGini) <- "giniImportance"
    machfile = paste(getwd(), "/randomForest/RF_giniImportance.csv", sep = "")
    write.csv(machGini, machfile, row.names = TRUE)
    
    machpvalue <- statRF$pimpTest$pvalue
    machfilep = paste(getwd(), "/randomForest/RF_pvalueImportance.csv", sep = "")
    write.csv(machpvalue, machfilep, row.names = TRUE)
    
    machfileMDS = paste(getwd(), "/randomForest/MultiDimensionalScalingPlot.pdf", sep = "")
    pdf(machfileMDS)
    mdsPlot <- mdsPlot(statRF$randomForest, statRF$pimpTest, Labels = Labels, slink = TRUE, slinkDat = getSlink)
    dev.off()
    machfileDim = paste(getwd(), "/randomForest/RF_mdsplotData.csv", sep = "")
    write.csv(mdsPlot, machfileDim, row.names = TRUE)
    
    cat("\n", "Permutation-based Gini importance measures")
    cat("\n", "The number of permutations: ", times)
    cat("\n", "Selected", nvarRF, "variables with top importance", "\n\n")
    
    machfileGini = paste(getwd(), "/randomForest/VarableImportancePlot.pdf", sep = "")
    pdf(machfileGini)
    giniPlot <- pvimPlot(statRF$randomForest, statRF$pimpTest, nvarRF = nvarRF)
    dev.off()
    selectedVar <- data.frame(giniImp = giniPlot$giniImp, pvalueImp = giniPlot$pvalueImp)
    selectedVar <- t(selectedVar)
    colnames(selectedVar) <- seq(1, nvarRF, 1)
    selectedVarPath = paste(getwd(), "/randomForest/RF_selectedVariablesData.csv", sep = "")
    write.csv(t(selectedVar), selectedVarPath, row.names = TRUE)
    print(selectedVar)
    
}



