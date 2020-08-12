outlier_stat <- function(pcx, pcy, scaling) {
    pwd.score = paste(getwd(), "/PCA_Data_", scaling, "/PCA_ScoreMatrix.csv", sep = "")
    Score <- read.csv(pwd.score, sep = ",", header = TRUE)
    Score.x <- Score[, 2:ncol(Score)]
    rownames(Score.x) <- Score[, 1]
    dx = scale(Score.x[, pcx], scale = FALSE)
    dy = scale(Score.x[, pcy], scale = FALSE)
    sumdxdx = sum(dx * dx)
    sumdydy = sum(dy * dy)
    sumdxdy = sum(dx * dy)
    theta = 0.5 * atan((2 * sumdxdy)/(sumdydy - sumdxdx))
    c = cos(theta)
    s = sin(theta)
    X = (c * dx) - (s * dy)
    Y = (s * dx) + (c * dy)
    varX = stats::var(X)
    varY = stats::var(Y)
    M = sqrt(varX)
    m = sqrt(varY)
    M95 = M * 3.03315
    m95 = m * 3.03315
    Fx = sqrt(abs((M95^2) - (m95^2)))
    Fy = Fx * tan(theta)
    F1 = c(-Fx, -Fy)
    F2 = c(Fx, Fy)
    one = matrix(rep(1, nrow(Score.x)), ncol = 1)
    F1.m = one %*% F1
    F2.m = one %*% F2
    # requireNamespace(pdist) library(pdist)
    Punti = cbind(dx, dy)
    dist1 = pdist::pdist(Punti[, 1:2], F1.m)
    dist2 = pdist::pdist(Punti[, 1:2], F2.m)
    dist1 = as.matrix(dist1)
    dist2 = as.matrix(dist2)
    D = matrix(dist1[, 1] + dist2[, 1], ncol = 1)
    v = M95 * 2
    outliers <- rownames(Score.x)[which(D > as.numeric(v))]
    n = length(outliers)
    cat("\n", "The following observations are calculated as outliers: ")
    if (n == 0) {
        print("No outliers were detected")
    } else {
        print(outliers)
    }
}
