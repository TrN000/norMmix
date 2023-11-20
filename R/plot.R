# plotting methods for norMmix objects


## colors:
nMmcols <- c(
    "#4363d8", "#f58231", "#800000", "#000075", "#ffe119",
    "#fabebe", "#e6beff", "#a9a9a9", "#ffffff", "#000000"
)
## chosen to be distinguishable and accessible for the colorblind,
## according to this site:
.colors_source <- "https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/"
## slightly rearranged, so that the first five colors stand out well
## on white background.

ellipsePts <- function(mu, sigma, npoints,
                       alpha = 0.05, r = sqrt(qchisq(1 - alpha, df = 2))) {
    stopifnot(is.matrix(sigma), length(mu) == 2L, dim(sigma) == 2L)
    theta <- seq(0, 2 * pi, length.out = npoints)
    es <- eigen(sigma)
    ## points
    rep(mu, each = npoints) -
        tcrossprod(
            r * cbind(cos(theta), sin(theta)), # v1
            es$vectors * sqrt(es$values)[c(1L, 1L, 2L, 2L)]
        ) # e1
}

## FIXME:  plot2d() <--> plotnd()  are *NOT* compatible in their defaults:
## =====
## 'npoints' should get the same *effective* default for p = 2
## 'border' had 'NA' for 2D, and polygon()'s default, 'NULL', for p > 2

ellipse_range <- function(nMm) {
    p <- nMm$dim
    k <- nMm$k

    ell_range <- matrix(c(Inf, -Inf), 2, p)

    # number of components
    for (i in 1:p) {
        for (j in 1:i) {
            if (i == j) next()
            ij <- c(i, j)
            pr <- norMmixProj(nMm, ij)
            for (comp in 1:k) {
                ell <- ellipsePts(pr$mu[, comp], pr$Sigma[, , comp], 100)
                range_i <- extendrange(ell[, 1], f = 0.25)
                range_j <- extendrange(ell[, 2], f = 0.25)
                if (ell_range[1, i] > range_i[1]) {
                    ell_range[1, i] <- range_i[1]
                }
                if (ell_range[2, i] < range_i[2]) {
                    ell_range[2, i] <- range_i[2]
                }
                if (ell_range[1, j] > range_j[1]) {
                    ell_range[1, j] <- range_j[1]
                }
                if (ell_range[2, j] < range_j[2]) {
                    ell_range[2, j] <- range_j[2]
                }
            }
        }
    }
    ell_range
}


## not exported; called from plot.norMmix() in 2D case, i.e. dim = 2
plot2d <- function(nMm, data = NULL,
                   plot = FALSE,
                   main = NULL,
                   sub = NULL,
                   type = "l", lty = 2, lwd = if (!is.null(data)) 2 else 1,
                   xlim = NULL, ylim = NULL, f.lim = 0.05,
                   npoints = 250, lab = FALSE,
                   col = nMmcols[1],
                   col.data = adjustcolor(par("col"), 1 / 2),
                   cex.data = par("cex"), pch.data = par("pch"),
                   fill = TRUE, fillcolor = col, border = NA,
                   ...) {
    w <- nMm$weight
    mu <- nMm$mu
    sig <- nMm$Sigma
    k <- nMm$k

    ## calculate smart values for xlim, ylim. if lims are given, variable xy still returned invisibly
    xy <- matrix(NA, k * npoints, 2)
    for (i in 1:k) {
        xy[(i - 1) * npoints + 1:npoints, ] <-
            ellipsePts(mu = mu[, i], sigma = sig[, , i], npoints = npoints)
    }
    if (is.null(xlim)) xlim <- extendrange(xy[, 1], f = f.lim)
    if (is.null(ylim)) ylim <- extendrange(xy[, 2], f = f.lim)

    if (fill) { ## determine fill color -- FIXME: use *correctly* below: different for components !!
        fco <- sapply(w, function(wj) adjustcolor(fillcolor, wj * 0.8 + 0.1))
    }

    ## add ellipses
    for (i in 1:k) {
        x <- ellipsePts(mu = mu[, i], sigma = sig[, , i], npoints = npoints)
        # if (i == 1) {
        if (plot && i == 1) {
            plot.default(x,
                type = type, xlim = xlim, ylim = ylim,
                main = main, sub = sub, lty = lty, col = col, ...
            )
        } else { # i > 1
            lines(x, type = type, lty = lty, col = col, ...)
        }

        if (fill) polygon(x, col = fco, border = border)
    }

    ## label components
    if (lab) {
        text(mu[1, ], mu[2, ], sprintf("comp %s", 1:k), adj = c(0.5, -4))
    }
    if (!is.null(data)) {
        points(data[, 1:2], cex = cex.data, col = col.data, pch = pch.data)
    }

    invisible(xy)
}


plotnd <- function(nMm, data = NULL,
                   main = NULL,
                   diag.panel = NULL,
                   ## for now, important arguments to mult.fig()
                   marP = rep(0, 4),
                   mgp = c(if (par("las") != 0) 2 else 1.5, 0.6, 0),
                   mar = marP + c(1.5, 2, .25, 0.5),
                   ...) {
    stopifnot(length(p <- nMm$dim) == 1, p >= 0)
    if (!is.null(diag.panel)) stopifnot(is.function(diag.panel), length(formals(diag.panel)) >= 1)

    xx <- ellipse_range(nMm)

    # off-diagonal panel function
    ## closure that captures nMm. we pass indices as 'x', 'y' values.
    ## determine the correct 2d projection per function call.
    simplot2d <- function(x, y, ij, data = NULL, main, sub, xlab, ylab, ...) { # swallow 'main', 'sub',...
        nm <- norMmixProj(nMm, ij)
        plot2d(nm, data = data, main = NULL, sub = NULL, xlab = NA, ylab = NA, plot = FALSE, ...)
    }

    pairs.indexed(
        xx,
        # maybe labels?,
        panel = simplot2d,
        diag.panel = diag.panel,
        gap = 0,
        # xlim = c(-3.0, 25.0), TODO: don't know how to set limits properly.
        # ylim = c(-3.0, 25.0),
    )
}


##' @title norMmix Projection on Coordinate Axes
##' @param nm norMmix obj
##' @param i # coordinate indices on which to project; integer vector, 1 <= i[k] <= p:= nm$dim
##' @return norMmix projected with reduced dimensions according to  ij
norMmixProj <- function(nm, ij) {
    nm$mu <- nm$ mu[ij, , drop = FALSE]
    nm$Sigma <- nm$ Sigma[ij, ij, , drop = FALSE]
    nm$dim <- length(ij)
    nm
}


plot.norMmixMLE <- function(x, y = NULL,
                            show.x = TRUE,
                            main = sprintf(
                                "norMmixMLE(*, model=\"%s\") fit to n=%d observations in %d dim.",
                                nm$model, x$nobs, nm$dim
                            ),
                            sub = paste0(
                                sprintf("log likelihood: %g; npar=%d", x$logLik, x$npar),
                                if (!is.null(opt <- x$optr)) {
                                    paste("; optim() counts:", named2char(opt$counts))
                                }
                            ),
                            cex.data = par("cex") / 4, pch.data = 4,
                            ...) {
    nm <- x$norMmix
    plot.norMmix(nm,
        data = if (show.x) x$x, # else NULL
        main = main, sub = sub, cex.data = cex.data, pch.data = pch.data, ...
    )
}


# plot function for norMmix objects
# \code{plot.norMmix} returns invisibly coordinates of bounding ellipses of distribution
plot.norMmix <- function(x, y = NULL, ...) {
    stopifnot(is.list(x), length(p <- x$dim) == 1)
    if (p == 2) {
        plot2d(x, ...)
    } else { ## if (p>2)
        plotnd(x, ...)
    } # replace with pairs function
}


### FIXME: Must document this ;  examples for *both* cases, ....
plot.fittednorMmix <- function(x, main = "unnamed", plotbest = FALSE, ...) {
    stopifnot(inherits(x, "fittednorMmix"))
    models <- x$models
    ## k <- x$k
    ## n <- x$n
    ## p <- x$p
    Bx <- BIC(x)
    bicmat <- Bx[[1]]
    best <- Bx[[2]]

    ### FIXME: should be able to plot *both* in one call ==> change UI (?!)
    if (!plotbest) {
        cl <- rainbow(length(models)) ## << FIXME! -- should be argument with *better* default

        matplot(bicmat, type = "l", xlab = "components", ylab = "BIC", col = cl, lty = 1:10, ...)
        legend("topright", models, fill = cl, lty = 1:10)
    } else {
        # TODO: like massplot (MM ??? )
        bk <- as.integer(best[1])
        bmodel <- best[2]
        plot(x$nMm[bk, bmodel][[1]]$norMmix, ...)
        points(x$x)
    }
    ## in both cases
    title(main = main)
    mtext(paste("best fit = ", best[1], best[2]))
}


pairs.indexed <- function(
        x, labels, panel = points, ..., horInd = 1:nc, verInd = 1:nc,
        lower.panel = panel, upper.panel = panel, diag.panel = NULL,
        text.panel = textPanel, label.pos = 0.5 + has.diag / 3, line.main = 3,
        cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1,
        log = "", horOdd = !row1attop, verOdd = !row1attop) {
    if (doText <- missing(text.panel) || is.function(text.panel)) {
        textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) {
            text(x, y, txt, cex = cex, font = font)
        }
    }
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, oma, ...) {
        xpd <- NA
        if (side %% 2L == 1L && xl[j]) {
            xpd <- FALSE
        }
        if (side %% 2L == 0L && yl[i]) {
            xpd <- FALSE
        }
        if (side %% 2L == 1L) {
            Axis(x, side = side, xpd = xpd, ...)
        } else {
            Axis(y, side = side, xpd = xpd, ...)
        }
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]])) {
                x[[i]] <- as.numeric(x[[i]])
            }
            if (!is.numeric(unclass(x[[i]]))) {
                stop("non-numeric argument to 'pairs'")
            }
        }
    } else if (!is.numeric(x)) {
        stop("non-numeric argument to 'pairs'")
    }
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) {
        lower.panel <- match.fun(lower.panel)
    }
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) {
        upper.panel <- match.fun(upper.panel)
    }
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) {
        diag.panel <- match.fun(diag.panel)
    }
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2L) {
        stop("only one column in the argument to 'pairs'")
    }
    if (!all(1L <= horInd & horInd <= nc)) {
        stop("invalid argument 'horInd'")
    }
    if (!all(1L <= verInd & verInd <= nc)) {
        stop("invalid argument 'verInd'")
    }
    if (doText) {
        if (missing(labels)) {
            labels <- colnames(x)
            if (is.null(labels)) {
                labels <- paste("var", 1L:nc)
            }
        } else if (is.null(labels)) {
            doText <- FALSE
        }
    }
    oma <- if ("oma" %in% nmdots) {
        dots$oma
    }
    main <- if ("main" %in% nmdots) {
        dots$main
    }
    if (is.null(oma)) {
        oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
    }
    opar <- par(mfcol = c(length(horInd), length(verInd)), mar = rep.int(
        gap / 2,
        4
    ), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    xl <- yl <- logical(nc)
    if (is.numeric(log)) {
        xl[log] <- yl[log] <- TRUE
    } else {
        xl[] <- grepl("x", log)
        yl[] <- grepl("y", log)
    }
    ni <- length(iSet <- if (row1attop) horInd else rev(horInd))
    nj <- length(jSet <- verInd)
    for (j in jSet) {
        for (i in iSet) {
            ij <- c(i, j)
            l <- paste0(if (xl[j]) {
                "x"
            } else {
                ""
            }, if (yl[i]) {
                "y"
            } else {
                ""
            })
            localPlot(x[, j], x[, i],
                xlab = "", ylab = "", axes = FALSE,
                type = "n", ..., log = l
            )
            if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
                box()
                j.odd <- (match(j, jSet) + horOdd) %% 2L
                i.odd <- (match(i, iSet) + verOdd) %% 2L
                if (i == iSet[1L] && (!j.odd || !has.upper || !has.lower)) {
                    localAxis(3L, x[, j], x[, i], ...)
                }
                if (i == iSet[ni] && (j.odd || !has.upper || !has.lower)) {
                    localAxis(1L, x[, j], x[, i], ...)
                }
                if (j == jSet[1L] && (!i.odd || !has.upper || !has.lower)) {
                    localAxis(2L, x[, j], x[, i], ...)
                }
                if (j == jSet[nj] && (i.odd || !has.upper || !has.lower)) {
                    localAxis(4L, x[, j], x[, i], ...)
                }
                mfg <- par("mfg")
                if (i == j) {
                    if (has.diag) {
                        localDiagPanel(as.vector(x[, i]), ij = ij, ...)
                    }
                    if (doText) {
                        par(usr = c(0, 1, 0, 1))
                        if (is.null(cex.labels)) {
                            l.wid <- strwidth(labels, "user")
                            cex.labels <- max(0.8, min(2, 0.9 / max(l.wid)))
                        }
                        xlp <- if (xl[i]) {
                            10^0.5
                        } else {
                            0.5
                        }
                        ylp <- if (yl[j]) {
                            10^label.pos
                        } else {
                            label.pos
                        }
                        text.panel(xlp, ylp, labels[i],
                            cex = cex.labels,
                            font = font.labels
                        )
                    }
                } else if (i < j) {
                    localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ij = ij, ...)
                } else {
                    localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), ij = ij, ...)
                }
                if (any(par("mfg") != mfg)) {
                    stop("the 'panel' function made a new plot")
                }
            } else {
                par(new = FALSE)
            }
        }
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) {
            dots$font.main
        } else {
            par("font.main")
        }
        cex.main <- if ("cex.main" %in% nmdots) {
            dots$cex.main
        } else {
            par("cex.main")
        }
        mtext(main, 3, line.main,
            outer = TRUE, at = 0.5, cex = cex.main,
            font = font.main
        )
    }
    invisible(NULL)
}
