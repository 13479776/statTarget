#' @name transX
#' @title transX for statTarget inputs
#' @description transX is to generate statTarget input file formats from 
#' Mass Spectrometry Data softwares, such as XCMS, MZmine2,SIEVE and SKYLINE.
#' @param data A transX objects. The output file from Mass Spectrometry Data 
#' softwares, such as '*.tsv' file from diffreport in XCMS software, 
#' '*.csv' file from SIEVE software,'.csv' file from SKYLINE software, 
#' '*.csv' file (Export to metaboAnalyst file) from MZmine2 software.
#' 
#' @param type The output file formats from Mass Spectrometry Data software, 
#' including 'XCMS' or 'xcms','MZmine2' or 'mzmine2','SIEVE' or 'sieve' and 
#' 'skyline' or 'SKYLINE'; Read-only .tsv file from diffreport in XCMS software
#' @return An output directory named 'statTargetDirectory'
#' @examples
#' datpath <- system.file('extdata',package = 'statTarget')
#' dataXcms <- paste(datpath,'xcmsOutput.tsv', sep='/')
#' dataSkyline <- paste(datpath,'skylineDemo.csv', sep='/')
#' transX(dataXcms,'xcms')
## transX(dataSkyline,'skyline')
#' @keywords XCMS MZmine2 SIEVE SKYLINE
#' @keywords inputs
#' @author Hemi Luan, hemi.luan@gmail.com 
#' @export 
transX <- function(data, type) {
    
    dirout.uni = paste(getwd(), "/statTargetDirectory/", sep = "")
    if (!file.exists("statTargetDirectory")) {
        dir.create(dirout.uni)
        # setwd(dirout.uni)
    }
    # dirout.w = paste(getwd(), '/statTargetDirectory/',type, sep='') dir.create(dirout.w)
    transX <- transCode(data = data, type = type)
    write.csv(transX$PhenoFile, paste(dirout.uni, "/metaFile_", type, ".csv", sep = ""), row.names = FALSE)
    write.csv(transX$ProfileFile, paste(dirout.uni, "/ProfileFile_", type, ".csv", sep = ""), row.names = FALSE)
    write.csv(transX$StatFile, paste(dirout.uni, "/StatFile_", type, ".csv", sep = ""), row.names = FALSE)
    write.csv(transX$info, paste(dirout.uni, "/SkyProtein", type, ".csv", sep = ""), row.names = FALSE)

    # setwd(dirout.uni)
    cat("Note: The input files have been generated for", type, ". Filling the missing info. please!\n")
}


transCode <- function(data, type) {
    # write.table(datR,'xcmsOutput_true.tsv',sep = '\t')
    if (type == "XCMS" | type == "xcms") {
        if (!grepl(".tsv", data) == TRUE) {
            stop("Read-only .tsv file from diffreport in XCMS software")
        }
        datR <- utils::read.delim(data)
        # ProfileFile
        ProfileFile <- cbind(datR$name, datR[, 15:ncol(datR)])
        colnames(ProfileFile)[1] <- c("name")
        sampleT <- c(colnames(ProfileFile[, 2:ncol(ProfileFile)]))
        ProteinID <- NA
        # PhenoFile
        PhenoFile <- as.data.frame(matrix(data = NA, nrow = length(sampleT), ncol = 4))
        colnames(PhenoFile) <- c("sample", "batch", "class", "order")
        PhenoFile$sample <- sampleT
        # statFile
        statF <- t(ProfileFile)
        colnames(statF) <- statF[1, ]
        statF <- data.frame(statF)
        statF <- cbind(PhenoFile$sample, PhenoFile$class, statF[-1, ])
        colnames(statF)[1:2] <- c("name", "group")
    }
    
    if (type == "SIEVE" | type == "sieve") {
        if (!grepl(".csv", data) == TRUE) {
            stop("Read-only .csv file from SIEVE software")
        }
        datR <- utils::read.csv(data, header = TRUE, sep = ",")
        # ProfileFile
        ProfileFile <- cbind(datR$CompID, datR[, 5:ncol(datR)])
        colnames(ProfileFile)[1] <- c("name")
        sampleT <- c(colnames(ProfileFile[, 2:ncol(ProfileFile)]))
        ProteinID <- NA
        # PhenoFile
        PhenoFile <- as.data.frame(matrix(data = NA, nrow = length(sampleT), ncol = 4))
        colnames(PhenoFile) <- c("sample", "batch", "class", "order")
        PhenoFile$sample <- sampleT
        # statFile
        statF <- t(ProfileFile)
        colnames(statF) <- statF[1, ]
        statF <- data.frame(statF)
        statF <- cbind(PhenoFile$sample, PhenoFile$class, statF[-1, ])
        colnames(statF)[1:2] <- c("name", "group")
    }
    
    if (type == "SKYLINE" | type == "skyline") {
        if (!grepl(".csv", data) == TRUE) {
            stop("Read-only .csv file from SKYLINE software")
        }
        datR <- utils::read.csv(data, header = TRUE, sep = ",")
        skyline <- function(x) {
            data <- as.data.frame(x)
            data$Area <- as.numeric(as.character(data$Area))
            data[is.na(data)] <- 0
            cols <- c("Precursor.Mz", "Product.Mz")
            datanew <- apply(data[, cols], 1, paste, collapse = "_")
            temp <- data[, c("Protein.Name", "Peptide.Sequence", "Replicate.Name", "Area")]
            datafile <- cbind(datanew, temp)
            colnames(datafile)[1] <- "name"
            uniID <- datafile[!duplicated(datafile[, c("name")]), ]
            cat(paste("Found", nrow(uniID), "targeted transitions"), "\n")
            resdat_temp <- plyr::daply(datafile, .(datanew, Replicate.Name), function(x) x$Area)
            resdat <- cbind(rownames(resdat_temp), resdat_temp)
            colnames(resdat)[1] <- "name"
            resdatjoin <- plyr::join(uniID, as.data.frame(resdat), by = "name")
            cat(paste("Conversion of skyline file was done !"), "\n")
            return(list(uniID = resdatjoin[, 1:3], area = resdatjoin[, -c(2, 3, 4, 5)]))
        }
        Skyoutput <- skyline(datR)
        # ProfileFile
        ProfileFile <- Skyoutput$area
        ProteinID <- Skyoutput$uniID
        sampleT <- c(colnames(ProfileFile[, 2:ncol(ProfileFile)]))
        # PhenoFile
        PhenoFile <- as.data.frame(matrix(data = NA, nrow = length(sampleT), ncol = 4))
        colnames(PhenoFile) <- c("sample", "batch", "class", "order")
        PhenoFile$sample <- sampleT
        # statFile
        statF <- t(ProfileFile)
        colnames(statF) <- statF[1, ]
        statF <- data.frame(statF)
        statF <- cbind(PhenoFile$sample, PhenoFile$class, statF[-1, ])
        colnames(statF)[1:2] <- c("name", "group")
        
        
        
        
    }
  
  if (type == "MZmine2" | type == "mzmine2") {
    if (!grepl(".csv", data) == TRUE) {
      stop("Read-only .csv file (Export to metaboAnalyst file) from MZmine2 software")
    }
    datR <- utils::read.csv(data)
    # ProfileFile
    ProfileFile <- as.data.frame(datR[-1, ])
    colnames(ProfileFile)[1] <- c("name")
    sampleT <- c(colnames(ProfileFile)[-1])
    ProteinID <- NA
    # PhenoFile
    PhenoFile <- as.data.frame(matrix(data = NA, nrow = c(length(ProfileFile) -1), ncol = 4))
    colnames(PhenoFile) <- c("sample", "batch", "class", "order")
    PhenoFile$sample <- sampleT
    # statFile
    statF <- t(ProfileFile)
    colnames(statF) <- statF[1, ]
    statF <- data.frame(statF)
    statF <- cbind(PhenoFile$sample, PhenoFile$class, statF[-1, ])
    colnames(statF)[1:2] <- c("name", "group")
  }
  
  dataOutput <- list(PhenoFile = PhenoFile, ProfileFile = ProfileFile, StatFile = statF, info = ProteinID)
    return(dataOutput)
}






