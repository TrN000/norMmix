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
                ell <- ellipsePts(pr$mu[, comp], pr$Sigma[, , comp], 1000)
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
                   add = FALSE,
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
        if (!add && i == 1) {
            plot.default(x,
                type = type, xlim = xlim, ylim = ylim,
                main = main, sub = sub, lty = lty, col = col, ...
            )
        } else { # add || i > 1
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
                   ...) {
    stopifnot(length(p <- nMm$dim) == 1, p >= 0)
    if (!is.null(diag.panel)) stopifnot(is.function(diag.panel), length(formals(diag.panel)) >= 1)

    xx <- as.jmatrix(ellipse_range(nMm))

    getJ <- function(.) attr(., "j")

    ## determine the correct 2d projection per function call through jmatrix trick.
    panel_2d <- function(x, y, data = NULL, ...) {
        i <- getJ(x)
        j <- getJ(y)
        nm <- norMmixProj(nMm, c(i, j))
        plot2d(nm, data = data, add=TRUE, ...)
    }

    pairs(
        xx,
        panel = panel_2d,
        diag.panel = diag.panel,
        gap = 0, # TODO: maybe more dynamic choice here.
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



as.jmatrix <- function(x) {
    stopifnot(length(d <- dim(x)) == 2)
    class(x) <- c("jmatrix", class(x))
    x
}

"[.jmatrix" <- function(x, i,j, drop=FALSE) {
    ret <- NextMethod()
    if(missing(i) && length(j) == 1L)
        structure(ret, j = j, class = "jvector")
    else
        structure(ret, class = class(x))
}

as.vector.jvector <- function(x, ...) x # no-op; needed as pairs.default calls  as.vector(x[,j])


