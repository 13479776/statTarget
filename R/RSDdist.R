RSDdist <- function(sample.rsd, sample.nor.rsd, QC.rsd, QC.nor.rsd) {
    # browser()
    colour1 <- NULL
    colour2 <- NULL
    colour1[(sample.nor.rsd/sample.rsd) > 1] <- "#D55E00"
    colour1[(sample.nor.rsd/sample.rsd) == 1] <- "grey48"
    colour1[(sample.nor.rsd/sample.rsd) < 1] <- "#009E73"
    
    colour2[(QC.nor.rsd/QC.rsd) > 1] <- "#D55E00"
    colour2[(QC.nor.rsd/QC.rsd) == 1] <- "grey48"
    colour2[(QC.nor.rsd/QC.rsd) < 1] <- "#009E73"
    
    s.rsd.up <- sum(colour1 == "#D55E00") * 100/length(colour1)
    s.rsd.no <- sum(colour1 == "grey48") * 100/length(colour1)
    s.rsd.down <- sum(colour1 == "#009E73") * 100/length(colour1)
    
    q.rsd.up <- sum(colour2 == "#D55E00") * 100/length(colour2)
    q.rsd.no <- sum(colour2 == "grey48") * 100/length(colour2)
    q.rsd.down <- sum(colour2 == "#009E73") * 100/length(colour2)
    
    dirout.w = paste(getwd(), "/statTarget/shiftCor/RSDresult", sep = "")
    
    pdf(file.path(dirout.w, "RSD variation.pdf"), width = 10, height = 7)
    layout(matrix(c(1, 2), ncol = 2))
    # par(mar=c(5,5,4,2))
    plot(sample.rsd, sample.nor.rsd, xlab = "RSD (Before correction)", ylab = "RSD (After correcion)", 
        col = colour1, cex.lab = 0.8, cex.axis = 0.8, main = "RSD variation of samples", cex.main = 1, 
        pch = 19)
    abline(0, 1, lwd = 1, lty = 2)
    abline(h = 30, lwd = 1, lty = 2)
    abline(v = 30, lwd = 1, lty = 2)
    
    legend("topleft", c(paste("Increased RSD after correction:", round(s.rsd.up, 2), "%"), paste("No changed RSD after correction:", 
        round(s.rsd.no, 2), "%"), paste("Decreased RSD after correction:", round(s.rsd.down, 2), "%")), 
        col = c("#D55E00", "grey48", "#009E73"), pch = 19, cex = 0.6)
    
    
    plot(QC.rsd, QC.nor.rsd, xlab = "RSD (Before correction)", ylab = "RSD (After correction)", col = colour2, 
        cex.lab = 0.8, cex.axis = 0.8, main = "RSD variation QC", cex.main = 1, pch = 19)
    
    abline(0, 1, lwd = 1, lty = 2)
    abline(h = 30, lwd = 1, lty = 2)
    abline(v = 30, lwd = 1, lty = 2)
    
    
    legend("topleft", c(paste("Increased RSDafter correction:", round(q.rsd.up, 2), "%"), paste("Not changed RSD after correction:", 
        round(q.rsd.no, 2), "%"), paste("Decreased RSD after correction:", round(q.rsd.down, 2), "%")), 
        col = c("#D55E00", "grey48", "#009E73"), pch = 19, cex = 0.6)
    dev.off()
    ## 
    s.rsd.dis <- sapply(seq(0, 1.9, 0.1), function(x) {
        sum(sample.rsd > x & sample.rsd <= x + 0.1)
    }) * 100/length(sample.rsd)
    s.nor.rsd.dis <- sapply(seq(0, 1.9, 0.1), function(x) {
        sum(sample.nor.rsd > x & sample.nor.rsd <= x + 0.1)
    }) * 100/length(sample.nor.rsd)
    q.rsd.dis <- sapply(seq(0, 1.9, 0.1), function(x) {
        sum(QC.rsd > x & QC.rsd <= x + 0.1)
    }) * 100/length(QC.rsd)
    q.nor.rsd.dis <- sapply(seq(0, 1.9, 0.1), function(x) {
        sum(QC.nor.rsd > x & QC.nor.rsd <= x + 0.1)
    }) * 100/length(QC.nor.rsd)
    
    rsd.dis <- rbind(s.rsd.dis, s.nor.rsd.dis, q.rsd.dis, q.nor.rsd.dis)
    colnames(rsd.dis) <- paste(paste(seq(0, 190, 10), seq(10, 200, 10), sep = "-"), "%", sep = "")
    rownames(rsd.dis) <- c("sample", "sample.nor", "QC", "QC.nor")
    # write.csv(rsd.dis,file.path(path,'RSD distribution.csv'))
    
    # par(new=TRUE)
    pdf(file.path(dirout.w, "RSD distribution.pdf"), width = 8, height = 6)
    layout(matrix(c(1, 2), ncol = 2))
    # par(mar=c(3,6,4,5))
    barp <- barplot(rsd.dis[1:2, ], horiz = TRUE, beside = TRUE, col = c("#D55E00", "#009E73"), border = c("#D55E00", 
        "#009E73"), names.arg = paste(seq(0, 190, 10), seq(10, 200, 10), sep = "-"), ylab = "RSD (%)", 
        las = 2, cex.lab = 0.6, cex.main = 0.8, cex.axis = 0.6, cex.names = 0.5, main = "Samples")
    text(rsd.dis[1:2, ] + 1, barp, labels = round(rsd.dis[1:2, ], 0), xpd = TRUE, cex = 0.3, col = "grey48")
    legend("topright", c("Before Correction", "After Correction"), col = c("#D55E00", "#009E73"), 
        lty = 0.5, pch = 15, cex = 0.75, horiz = FALSE, bty = "n")
    # par(mar=c(5,6,4,5))
    barp <- barplot(rsd.dis[3:4, ], horiz = TRUE, beside = TRUE, col = c("#D55E00", "#009E73"), border = c("#D55E00", 
        "#009E73"), names.arg = paste(seq(0, 190, 10), seq(10, 200, 10), sep = "-"), ylab = "RSD (%)", 
        las = 2, cex.lab = 0.6, cex.main = 0.8, cex.axis = 0.6, cex.names = 0.5, main = "Quality Controls (QC)")
    text(rsd.dis[3:4, ] + 1, barp, labels = round(rsd.dis[3:4, ], 0), xpd = TRUE, cex = 0.3, col = "grey48")
    legend("topright", c("Before Correction", "After Correction"), col = c("#D55E00", "#009E73"), 
        lty = 0.5, pch = 15, cex = 0.75, horiz = FALSE, bty = "n")
    # dev.copy2pdf(file.path(dirout.w,'RSD distribution.pdf'),width=8,height=6)
    dev.off()
}
