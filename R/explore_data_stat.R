explore_data_stat <- function(file, scaling, scal = TRUE, normalize = FALSE, imputation = FALSE, imput) {
    comp = read.csv(file, sep = ",", header = TRUE)
    comp.x = comp[, 3:ncol(comp)]
    comp.x = cbind(comp[, 2], comp[, 1], comp.x)
    x <- comp.x
    x.x <- x[, 3:ncol(x)]
    rownames(x.x) <- x[, 2]
    dirout = paste(getwd(), "/scaleData_", scaling, "/", sep = "")
    dir.create(dirout)
    x.x = cbind(comp[, 2], x.x)
    write.csv(x.x, paste(dirout, "CorrectedTable.csv", sep = ""))
    k = matrix(x[, 1], ncol = 1)
    K = paste(getwd(), "/scaleData_", scaling, "/class.csv", sep = "")
    write.csv(k, K)
    # scaling
    if (scal) {
        if (scaling == "Pareto" | scaling == "pareto" | scaling == "P" | scaling == "p") {
            pwd.n = paste(getwd(), "/scaleData_", scaling, "/CorrectedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 3:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            # x.areanorm.tc[x.areanorm.tc == 0 ] <- runif(1, 0, 0.0000000001)
            
            all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
            all.sd[all.sd == 0] <- runif(1, 0, 1e-10)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), ncol = 1)
            all.sdm = uni.exp.all %*% all.sd
            all.sqsd = sqrt(all.sdm)
            all.pareto <- x.areanorm.tc/all.sqsd
            write.csv(all.pareto, paste(dirout, "/ProcessedTable.csv", sep = ""))
        } else if (scaling == "Auto" | scaling == "auto" | scaling == "A" | scaling == "a") {
            pwd.n = paste(getwd(), "/scaleData_", scaling, "/CorrectedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 3:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
            all.sd[all.sd == 0] <- runif(1, 0, 1e-10)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), ncol = 1)
            all.sdm = uni.exp.all %*% all.sd
            all.auto <- x.areanorm.tc/all.sdm
            write.csv(all.auto, paste(dirout, "/ProcessedTable.csv", sep = ""))
        } else if (scaling == "Vast" | scaling == "vast" | scaling == "V" | scaling == "v") {
            pwd.n = paste(getwd(), "/scaleData_", scaling, "/CorrectedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 3:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            all.sd <- matrix(apply(x.areanorm.tc, 2, sd), nrow = 1)
            all.sd[all.sd == 0] <- runif(1, 0, 1e-10)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), ncol = 1)
            all.sdm = uni.exp.all %*% all.sd
            sdm2 = all.sdm^2
            colm = matrix(colMeans(x.x), nrow = 1)
            colm.m = uni.exp.all %*% colm
            num = x.areanorm.tc * colm.m
            vast = num/sdm2
            write.csv(vast, paste(dirout, "/ProcessedTable.csv", sep = ""))
        } else if (scaling == "Range" | scaling == "range" | scaling == "R" | scaling == "r") {
            pwd.n = paste(getwd(), "/scaleData_", scaling, "/CorrectedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 3:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            range = c()
            for (i in 1:ncol(x.x)) {
                den = c()
                den = max(x.x[, i]) - min(x.x[, i])
                range = matrix(c(range, den), nrow = 1)
            }
            range[range == 0] <- runif(1, 0, 1e-10)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), ncol = 1)
            range.m = uni.exp.all %*% range
            all.range = x.areanorm.tc/range.m
            write.csv(all.range, paste(dirout, "/ProcessedTable.csv", sep = ""))
        } else if (scaling == "Median" | scaling == "median" | scaling == "M" | scaling == "m") {
            pwd.n = paste(getwd(), "/scaleData_", scaling, "/CorrectedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 3:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.areanorm.tc <- scale(x.x, center = TRUE, scale = FALSE)
            all.med <- matrix(apply(x.areanorm.tc, 2, median), nrow = 1)
            uni.exp.all = matrix(rep(1, nrow(x.areanorm.tc)), ncol = 1)
            all.sdm = uni.exp.all %*% all.med
            all.med <- x.areanorm.tc/all.sdm
            write.csv(all.med, paste(dirout, "/ProcessedTable.csv", sep = ""))
        } else if (scaling == "NONE" | scaling == "none" | scaling == "N" | scaling == "n") {
            pwd.n = paste(getwd(), "/scaleData_", scaling, "/CorrectedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 3:ncol(x)]
            rownames(x.x) <- x[, 1]
            all.med <- x.x
            write.csv(all.med, paste(dirout, "/ProcessedTable.csv", sep = ""))
            
        } else if (scaling == "Center" | scaling == "center" | scaling == "C" | scaling == "c") {
            pwd.n = paste(getwd(), "/scaleData_", scaling, "/CorrectedTable.csv", sep = "")
            x <- read.csv(pwd.n, sep = ",", header = TRUE)
            x.x <- x[, 3:ncol(x)]
            rownames(x.x) <- x[, 1]
            x.c = scale(x.x, scale = FALSE)
            write.csv(x.c, paste(dirout, "/ProcessedTable.csv", sep = ""))
        }
    }
    pwd.scal = paste(getwd(), "/scaleData_", scaling, "/ProcessedTable.csv", sep = "")
    x <- read.csv(pwd.scal, sep = ",", header = TRUE)
    x.x <- x[, 2:ncol(x)]
    rownames(x.x) <- x[, 1]
    pc.all <- prcomp(x.x, center = FALSE, scale = FALSE)
    cat("\n", "PCA Model Summary")
    cat("\n", nrow(comp), "samples x", ncol(comp) - 2, "variables")
    cat("\n", "Variance Explained of PCA Model: ", "\n\n")
    summPca <- summary(pc.all)$importance
    if (dim(summPca)[2] < 5) {
        print(summPca)
    } else {
        print(summPca[, 1:5])
    }
    p.v <- matrix(((pc.all$sdev^2)/(sum(pc.all$sdev^2))), ncol = 1)
    p.i <- round(p.v * 100, 1)
    p.z <- matrix(1, nrow(p.i), 1)
    p.f <- cbind(p.i, p.z)
    dirout.pca = paste(getwd(), "/PCA_Data_", scaling, "/", sep = "")
    dir.create(dirout.pca)
    write.csv(p.f, paste(dirout.pca, "PCA_P", sep = ""))
    write.csv(pc.all$x, paste(dirout.pca, "PCA_ScoreMatrix.csv", sep = ""))
    
    write.csv(t(summPca), paste(dirout.pca, "PCA_explainVariance.csv", sep = ""))
    
    write.csv(pc.all$rotation, paste(dirout.pca, "PCA_LoadingsMatrix.csv", sep = ""))
    pwd.score = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_ScoreMatrix.csv", sep = "")
    Score <- read.csv(pwd.score, sep = ",", header = TRUE)
    Score.x <- Score[, 2:ncol(Score)]
    rownames(Score.x) <- Score[, 1]
    pwd.load = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_LoadingsMatrix.csv", sep = "")
    Loading <- read.csv(pwd.load, sep = ",", header = TRUE)
    Loading.x <- Loading[, 2:ncol(Loading)]
    rownames(Loading.x) <- Loading[, 1]
    pwd.pvar = paste(getwd(), "/PCA_Data_", scaling, "/", "PCA_P", sep = "")
    Pvar <- read.csv(pwd.pvar, sep = ",", header = TRUE)
    Pvar.x <- Pvar[, 2:ncol(Pvar)]
    rownames(Pvar.x) <- Pvar[, 1]
    scree = paste(dirout.pca, "Screeplot", scaling, ".pdf", sep = "")
    pdf(scree)
    barplot(Pvar.x[, 1], xlab = "Principal Components", ylab = "Proportion of Variance explained", 
        main = "Screeplot", ylim = c(0, 100))
    dev.off()
    pairs = paste(dirout.pca, "First_3_Components_", scaling, ".pdf", sep = "")
    pdf(pairs)
    tutticolors = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, "rosybrown4", "green4", "navy", "purple2", "orange", 
        "pink", "chocolate2", "coral3", "khaki3", "thistle", "turquoise3", "palegreen1", "moccasin", 
        "olivedrab3", "azure4", "gold3", "deeppink"), ncol = 1)
    col = c()
    for (i in 1:nrow(k)) {
        col = c(col, tutticolors[k[i, ], ])
    }
    pairs = c()
    if (ncol(Score.x) >= 3) {
        pairs = c(3)
    } else {
        pairs = c(ncol(Score.x))
    }
    pairs(Score.x[, 1:pairs], col = col)
    dev.off()
}
