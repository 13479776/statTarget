dataEllipse_sT <- function(x, y, groups, group.labels = group.levels, ellipse.label, weights, log = "", 
    levels = c(0.5, 0.95), center.pch = 19, center.cex = 1.5, draw = TRUE, plot.points = draw, add = !plot.points, 
    segments = 51, robust = FALSE, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), col = if (missing(groups)) palette()[1:2] else palette()[1:length(group.levels)], 
    pch = if (missing(groups)) 1 else seq(group.levels), lwd = 2, fill = FALSE, fill.alpha = 0.3, 
    grid = TRUE, labels, id.method = "mahal", id.n = if (id.method[1] == "identify") Inf else 0, id.cex = 1, 
    id.col = if (missing(groups)) palette()[1] else palette()(1:length(groups)), ...) {
    label.ellipse <- function(ellipse, label, col, ...) {
        # Ellipses (orignally by J. Fox and G. Monette) This sub-function from Michael Friendly position
        # label above top right
        if (cor(ellipse)[1, 2] >= 0) {
            index <- which.max(ellipse[, 2])
            x <- ellipse[index, 1] + 0.5 * strwidth(label)
            y <- ellipse[index, 2] + 0.5 * strheight("A")
            adj <- c(1, 0)
        } else {
            # position label below bot left
            index <- which.min(ellipse[, 2])
            x <- ellipse[index, 1] - 0.5 * strwidth(label)
            y <- ellipse[index, 2] - 0.5 * strheight("A")
            adj <- c(0, 1)
        }
        text(x, y, label, adj = adj, col = col, ...)
    }
    if (missing(y)) {
        if (is.matrix(x) && ncol(x) == 2) {
            if (missing(xlab)) 
                xlab <- colnames(x)[1]
            if (missing(ylab)) 
                ylab <- colnames(x)[2]
            y <- x[, 2]
            x <- x[, 1]
        } else stop("x and y must be vectors, or x must be a 2 column matrix")
    } else if (!(is.vector(x) && is.vector(y) && length(x) == length(y))) 
        stop("x and y must be vectors of the same length")
    if (missing(weights)) 
        weights <- rep(1, length(x))
    if (length(weights) != length(x)) 
        stop("weights must be of the same length as x and y")
    if (!missing(groups)) {
        xlab
        ylab
        if (!is.factor(groups)) 
            stop("groups must be a factor")
        if (!(length(groups) == length(x))) 
            stop("groups, x, and y must all be of the same length")
        if (missing(labels)) 
            labels <- seq(length(x))
        valid <- complete.cases(x, y, groups)
        x <- x[valid]
        y <- y[valid]
        weights <- weights[valid]
        groups <- groups[valid]
        labels <- labels[valid]
        group.levels <- levels(groups)
        col <- col[!is.na(col)]
        if (length(col) < length(group.levels)) 
            stop("too few colors for number of groups")
        result <- vector(length(group.levels), mode = "list")
        names(result) <- group.levels
        if (draw) {
            if (!add) {
                plot(x, y, type = "n", xlab = xlab, ylab = ylab, ...)
                if (grid) {
                  grid(lty = 1, equilogs = FALSE)
                  box()
                }
            }
        }
        for (lev in 1:length(group.levels)) {
            level <- group.levels[lev]
            sel <- groups == level
            result[[lev]] <- dataEllipse_sT(x[sel], y[sel], weights = weights[sel], log = log, levels = levels, 
                center.pch = center.pch, center.cex = center.cex, draw = draw, plot.points = plot.points, 
                add = TRUE, segments = segments, robust = robust, col = rep(col[lev], 2), pch = pch[lev], 
                lwd = lwd, fill = fill, fill.alpha = fill.alpha, labels = labels[sel], id.method = id.method, 
                id.n = id.n, id.cex = id.cex, id.col = col[lev], ellipse.label = group.labels[lev], 
                ...)
        }
        return(invisible(result))
    }
    if (length(col) == 1) 
        col <- rep(col, 2)
    if (draw) {
        if (!add) {
            plot(x, y, type = "n", xlab = xlab, ylab = ylab, ...)
            if (grid) {
                grid(lty = 1, equilogs = FALSE)
                box()
            }
        }
        if (plot.points) 
            points(x, y, col = col[1], pch = pch[1], ...)
    }
    dfn <- 2
    dfd <- length(x) - 1
    if (robust) {
        use <- weights > 0
        v <- cov_trob_sT(cbind(x[use], y[use]), wt = weights[use])
        shape <- v$cov
        center <- v$center
    } else {
        v <- cov.wt(cbind(x, y), wt = weights)
        shape <- v$cov
        center <- v$center
    }
    result <- vector("list", length = length(levels))
    names(result) <- levels
    for (i in seq(along = levels)) {
        level <- levels[i]
        radius <- sqrt(dfn * qf(level, dfn, dfd))
        result[[i]] <- ellipse_sT(center, shape, radius, log = log, center.pch = center.pch, center.cex = center.cex, 
            segments = segments, col = col[2], lwd = lwd, fill = fill, fill.alpha = fill.alpha, draw = draw, 
            ...)
        if (!missing(ellipse.label)) {
            lab <- if (length(ellipse.label) < i) 
                ellipse.label[1] else ellipse.label[i]
            label.ellipse(result[[i]], lab, col[2], ...)
        }
    }
    if (missing(labels)) 
        labels <- seq(length(x))
    if (draw) 
        showLabels_sT(x, y, labels = labels, id.method = id.method, id.n = id.n, id.cex = id.cex, 
            id.col = id.col)
    invisible(if (length(levels) == 1) result[[1]] else result)
}
