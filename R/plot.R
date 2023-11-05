# plotting methods for norMmix objects


## colors:
nMmcols <- c("#4363d8", "#f58231", "#800000", "#000075", "#ffe119",
             "#fabebe", "#e6beff", "#a9a9a9", "#ffffff", "#000000")
## chosen to be distinguishable and accessible for the colorblind,
## according to this site:
## https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
## slightly rearranged, so that the first five colors stand out well
## on white background.

## MM:  the part of mixtools::ellipse() we need is trivial:
ellipsePts <- function(mu, sigma, npoints,
                       alpha = 0.05, r = sqrt(qchisq(1 - alpha, df = 2))) {
    stopifnot(is.matrix(sigma), length(mu) == 2L, dim(sigma) == 2L)
    theta <- seq(0, 2 * pi, length.out = npoints)
    es <- eigen(sigma)
    ## points
    rep(mu, each = npoints) -
        tcrossprod(r * cbind(cos(theta), sin(theta)), # v1
                   es$vectors * sqrt(es$values)[c(1L, 1L, 2L, 2L)]) # e1
}

## FIXME:  plot2d() <--> plotnd()  are *NOT* compatible in their defaults:
## =====
## 'npoints' should get the same *effective* default for p = 2
## 'border' had 'NA' for 2D, and polygon()'s default, 'NULL', for p > 2


## not exported; called from plot.norMmix() in 2D case, i.e. dim = 2
plot2d <- function(nMm, data = NULL,
                   main = deparse(sys.call()), # <- rather main = NULL ?
                   sub = NULL,
                   type = "l", lty = 2, lwd  =  if (!is.null(data)) 2 else 1,
                   xlim = NULL, ylim = NULL, f.lim = 0.05,
                   npoints = 250, lab = FALSE,
                   col = nMmcols[1],
                   col.data  =  adjustcolor(par("col"), 1 / 2),
                   cex.data  =  par("cex"), pch.data  =  par("pch"),
                   fill = TRUE, fillcolor = col, border = NA,
                   ...)  {
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

    if (fill) ## determine fill color -- FIXME: use *correctly* below: different for components !!
        fco <- sapply(w, function(wj) adjustcolor(fillcolor, wj * 0.8 + 0.1))

    ## add ellipses
    for (i in 1:k) {
        x <- ellipsePts(mu = mu[, i], sigma = sig[, , i], npoints = npoints)
        if (i == 1)
            plot.default(x, type = type, xlim = xlim, ylim = ylim,
                         main = main, sub = sub, lty = lty, col = col, ...)
        else # i > 1
            lines(x, type = type, lty = lty, col = col, ...)

        if (fill) polygon(x, col = fco, border = border)
    }

    ## label components
    if (lab)
        text(mu[1, ], mu[2, ], sprintf("comp %s", 1:k), adj = c(0.5, -4))
    if (!is.null(data))
        points(data[, 1:2], cex = cex.data, col = col.data, pch = pch.data)

    invisible(xy)
}


plotnd <- function(nMm, data = NULL,
                   ###
                   ### FIXME: (see ../TODO.md) use pairs() with correct panel = function(.).....
                   ###
                   main = NULL, ## better than main=deparse(sys.call()),
                   diag.panel = NULL,
                   ## for now, important arguments to mult.fig()
                   marP = rep(0, 4),
                   mgp = c(if (par("las") != 0) 2 else 1.5, 0.6, 0),
                   mar = marP + c(1.5, 2, .25, 0.5),
                   ...) {
    stopifnot(length(p <- nMm$dim) == 1, p >= 0)

    # preallocating result (a "matrix" list):
    pts <- vector("list", p * p)
    dim(pts) <- c(p, p)
    if (p == 0) return(invisible(pts))


    # defining diagonal panel
    if (is.null(diag.panel)) {
        diag.panel <- function(i, p) {
            frame()
            box()
            text(0.5, 0.5, paste("var", i, cex = p)) # rather cex = 2 (??)
        }
    } else {# a function with at least one argument :
        stopifnot(is.function(diag.panel), length(formals(diag.panel)) >= 2)
    }

    simplot2d <- function(nm, data, main, sub, xlab, ylab, ...) { # swallow 'main', 'sub',...
        plot2d(nm, data = data, main = NULL, sub = NULL, xlab = NA, ylab = NA, ...)
    }

    # setting up the plot _______ FIXME Use  pairs() with panel in the future !!
    opar <- mult.fig(mfcol = c(p, p), marP = marP, mgp = mgp, mar = mar, main = main) $ opar
    on.exit(par(opar))
    for (i in 1:p) {
        for (j in 1:p) {
            if (j == i) { # diagonal panel
                diag.panel(i, p)
            } else { # off diag panels, plot 2-dim projections using ellipsePts()
                pts[i, j][[1L]] <- simplot2d(norMmixProj(nMm, c(i, j)), data = data[, c(i, j)], ...)
            }
        }
    }
    invisible(pts)
}


##' @title norMmix Projection on Coordinate Axes
##' @param nm norMmix obj
##' @param i # coordinate indices on which to project; integer vector, 1 <= i[k] <= p:= nm$dim
##' @return norMmix projected with reduced dimensions according to  ij
norMmixProj <- function(nm, ij) {
    nm$mu    <-  nm$ mu   [ij,     , drop = FALSE]
    nm$Sigma <-  nm$ Sigma[ij, ij, , drop = FALSE]
    nm$dim <- length(ij)
    nm
}


plot.norMmixMLE <- function(x, y = NULL,
                            show.x = TRUE,
                            main = sprintf("norMmixMLE(*, model=\"%s\") fit to n=%d observations in %d dim.",
                                           nm$model, x$nobs, nm$dim),
                            sub = paste0(sprintf("log likelihood: %g; npar=%d", x$logLik, x$npar),
                                         if (!is.null(opt <- x$optr))
                                             paste("; optim() counts:", named2char(opt$counts))),
                            cex.data = par("cex") / 4, pch.data = 4,
                            ...) {
    nm <- x$norMmix
    plot.norMmix(nm, data = if (show.x) x$x, # else NULL
                 main = main, sub = sub, cex.data = cex.data, pch.data = pch.data, ...)
}


# plot function for norMmix objects
# \code{plot.norMmix} returns invisibly coordinates of bounding ellipses of distribution
plot.norMmix <- function(x, y = NULL, ...) {
    stopifnot(is.list(x), length(p <- x$dim) == 1)
    if (p == 2)
        plot2d(x, ...)
    else ## if (p>2)
        plotnd(x, ...)
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
    best   <- Bx[[2]]

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
