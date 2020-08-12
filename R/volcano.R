## FoldChange value is calculated from the raw data (dataSummary file).
volcano <- function(file, upper.lim, lower.lim, sig.lim) {
    pwdfile = paste(getwd(), "/Univariate/DataTable.csv", sep = "")
    file = pwdfile
    x <- read.csv(file, sep = ",", header = TRUE)
    x.x = x[, 3:ncol(x)]
    rownames(x.x) = x[, 2]
    k = matrix(x[, 1], ncol = 1)
    slink = paste(getwd(), "/DataPretreatment", "/slink.csv", sep = "")
    slink = read.csv(slink, header = TRUE)
    x.n = cbind(k, x.x)
    sorted = x.n[order(x.n[, 1]), ]
    sorted.x = as.matrix(sorted[, -1], ncol = ncol(sorted) - 1)
    g = c()
    for (i in 1:nrow(sorted)) {
        if (any(g == sorted[i, 1])) {
            g = g
        } else {
            g = matrix(c(g, sorted[i, 1]), ncol = 1)
        }
    }
    NoF = nrow(g)
    dirout.fc = paste(getwd(), "/Univariate/Fold_Changes/", sep = "")
    dir.create(dirout.fc)
    dirout.vol = paste(getwd(), "/Univariate/Volcano_Plots/", sep = "")
    dir.create(dirout.vol)
    for (i in 1:NoF) {
        for (j in 1:NoF) {
            if (i < j) {
                ni = paste("dataSumm_", ExcName(i, slink), ".csv", sep = "")
                nj = paste("dataSumm_", ExcName(j, slink), ".csv", sep = "")
                # pwdi = paste(getwd(), '/Univariate/Groups/', ni, sep='') pwdj = paste(getwd(),
                # '/Univariate/Groups/', nj, sep='')
                pwdi = paste(getwd(), "/dataSummary/", ni, sep = "")
                pwdj = paste(getwd(), "/dataSummary/", nj, sep = "")
                pv = paste(getwd(), "/Univariate/Pvalues/Pvalues_", ExcName(i, slink), "vs", ExcName(j, 
                  slink), ".csv", sep = "")
                
                
                I = read.csv(pwdi, header = TRUE)
                J = read.csv(pwdj, header = TRUE)
                mI = I[, "Mean"]
                mJ = J[, "Mean"]
                
                FC = matrix(abs(mI/mJ), ncol = 1)
                rownames(FC) = I[, 1]
                colnames(FC) <- "foldChanges"
                fc.csvfile = paste("Fold_Change_", ExcName(i, slink), "-to-", ExcName(j, slink), ".csv", 
                  sep = "")
                write.csv(FC, paste(dirout.fc, fc.csvfile, sep = ""))
                PV = read.csv(pv, header = TRUE)
                if("pvalue" %in% colnames(PV)) adjCheck <- "pvalue" else adjCheck <- "adjust.pvalue"
                PV = matrix(PV[, -1], ncol = 1)
                PV[is.na(PV)] <- 1
                logfc = log2(FC)
                # logfc = glog(FC,2)
                logpv = -log10(PV)
                # logfc = FC
                
                upper <- log2(upper.lim)
                lower <- log2(lower.lim)
                # upper <- upper.lim lower <- lower.lim
                
                
                sig <- -log10(sig.lim)
                colpv = matrix(rep(NA, nrow(PV)), ncol = 1)
                rownames(colpv) <- rownames(logfc)
                for (p in 1:nrow(PV)) {
                  if (logfc[p, ] < lower & logpv[p, ] > sig) {
                    colpv[p, ] = "navy"
                  } else if (logfc[p, ] > upper & logpv[p, ] > sig) {
                    colpv[p, ] = "#D55E00"
                  } else {
                    colpv[p, ] = "dark grey"
                  }
                }
                
                max.fc = 1.3 * max((abs(logfc)))
                V = paste(dirout.vol, "VolcanoPlot_", ExcName(i, slink), "-to-", ExcName(j, slink), 
                  ".pdf", sep = "")
                pospv = matrix(rep(NA, nrow(PV)), ncol = 1)
                for (p in 1:nrow(PV)) {
                  if (logfc[p, ] < lower) {
                    pospv[p, ] = 2L
                  } else {
                    pospv[p, ] = 4L
                  }
                }
                
                pdf(V)
                graphics::plot(logfc, logpv, col = colpv, cex = 1.2, pch = 19, cex.sub = 0.8, xlim = c(-max.fc, 
                  max.fc), xlab = "log2 (Fold Change)", ylab = paste("-log10(",adjCheck,")",sep = ""), main = paste("Volcano Plot ", 
                  ExcName(i, slink), "-to-", ExcName(j, slink), sep = ""), sub = paste("(Blue dot, down-regulated; red dot, up-regulated; ", 
                  "Pvalue < ", sig.lim, ", FC > ", upper.lim, " or < ", lower.lim, ")", sep = ""))
                
                
                voldata <- cbind(logfc, logpv)
                
                if (length(colnames(sorted.x)[grep("navy", colpv)]) > 0 & length(colnames(sorted.x)[grep("#D55E00", 
                  colpv)]) > 0) {
                  text(voldata[grep("navy", colpv), 1], voldata[grep("navy", colpv), 2], labels = rownames(voldata)[grep("navy", 
                    colpv)], cex = 0.6, pos = pospv[grep("navy", colpv)], col = colpv[grep("navy", 
                    colpv)])
                  text(voldata[grep("#D55E00", colpv), 1], voldata[grep("#D55E00", colpv), 2], labels = rownames(voldata)[grep("#D55E00", 
                    colpv)], cex = 0.6, pos = pospv[grep("#D55E00", colpv)], col = colpv[grep("#D55E00", 
                    colpv)])
                } else if (length(colnames(sorted.x)[grep("navy", colpv)]) > 0 & length(colnames(sorted.x)[grep("#D55E00", 
                  colpv)]) == 0) {
                  text(voldata[grep("navy", colpv), 1], voldata[grep("navy", colpv), 2], labels = rownames(voldata)[grep("navy", 
                    colpv)], cex = 0.6, pos = pospv[grep("navy", colpv)], col = colpv[grep("navy", 
                    colpv)])
                } else if (length(colnames(sorted.x)[grep("navy", colpv)]) == 0 & length(colnames(sorted.x)[grep("#D55E00", 
                  colpv)]) > 0) {
                  text(voldata[grep("#D55E00", colpv), 1], voldata[grep("#D55E00", colpv), 2], labels = rownames(voldata)[grep("#D55E00", 
                    colpv)], cex = 0.6, pos = pospv[grep("#D55E00", colpv)], col = colpv[grep("#D55E00", 
                    colpv)])
                }
                
                
                
                if (length(colnames(sorted.x)[grep("navy", colpv)]) > 0 | length(colnames(sorted.x)[grep("#D55E00", 
                  colpv)]) > 0) {
                  
                  Vol_sig = paste(dirout.vol, "VolcanoMarker_", ExcName(i, slink), "-to-", ExcName(j, 
                    slink), ".csv", sep = "")
                  dirout.name = c(rownames(voldata)[grep("navy", colpv)], rownames(voldata)[grep("#D55E00", 
                    colpv)])
                  
                  dirout.name <- as.data.frame(dirout.name)
                  colnames(dirout.name) <- c("Diff_VolcanoMarker")
                  write.csv(dirout.name, Vol_sig, row.names = FALSE)
                } else {
                  NULL
                }
                axis(2, at = c(-1, 150), pos = c(lower, 0), col = "blue", lwd = 1, lty = 2)
                axis(2, at = c(-1, 150), pos = c(upper, 0), col = "blue", lwd = 1, lty = 2)
                axis(1, at = c(-150, 150), pos = c(sig, 0), col = "blue", lwd = 1, lty = 2)
                dev.off()
                
            }
        }
    }
}
