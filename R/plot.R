### plotting methods for norMmix objects


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
                       alpha = 0.05, r = sqrt(qchisq(1-alpha, df=2))) {
    stopifnot(is.matrix(sigma), length(mu) == 2L, dim(sigma) == 2L)
    theta <- seq(0, 2*pi, len=npoints)
    es <- eigen(sigma)
    ## points
    rep(mu, each=npoints) -
        tcrossprod(r * cbind(cos(theta), sin(theta)), # v1
                   es$vectors * sqrt(es$values)[c(1L,1L,2L,2L)]) # e1
}

## FIXME:  plot2d() <--> plotnd()  are *NOT* compatible in their defaults:
## =====
## 'npoints' should get the same *effective* default for p=2
## 'border' had 'NA' for 2D, and polygon()'s default, 'NULL', for p > 2


## not exported; called from plot.norMmix() in 2D case, i.e. dim=2
plot2d <- function(nMm, data=NULL , type="l", lty=2,
                   xlim=NULL, ylim=NULL, f.lim=0.05,
                   newWindow=TRUE, npoints=250, lab=FALSE,
                   col=nMmcols[1],
                   fill=TRUE, fillcolor=nMmcols[1], border=NA,
	           ...)  {
    w <- nMm$weight
    mu <- nMm$mu
    sig <- nMm$Sigma
    k <- nMm$k

    ## calculate smart values for xlim, ylim. if lims are given, variable xy still returned invisibly
    xy <- matrix(NA, k*npoints, 2)
    for (i in 1:k) {
        xy[(i-1)*npoints + 1:npoints, ] <-
            ellipsePts(mu=mu[,i], sigma=sig[,,i], npoints=npoints)
    }
    if (is.null(xlim)) xlim <- extendrange(xy[,1], f=f.lim)
    if (is.null(ylim)) ylim <- extendrange(xy[,2], f=f.lim)

    if(fill) ## determine fill color -- FIXME: use *correctly* below: different for components !!
        fco <- sapply(w, function(wj) adjustcolor(fillcolor, wj*0.8 + 0.1))

    ## add ellipses
    for (i in 1:k) {
        x <- ellipsePts(mu=mu[,i], sigma=sig[,,i], npoints=npoints)
        if(i == 1)
            plot(x, type=type, xlim=xlim, ylim=ylim, lty=lty, col=col, ...)
        else # i > 1
            points(x, type=type, lty=lty, col=col, ...)

        if (fill) polygon(x, col=fco, border=border)
    }

    ## label components
    if (lab)
        text(mu[1,], mu[2,], sprintf("comp %s", 1:k), adj=c(0.5,-4))
    if (!is.null(data))
        points(data[, 1:2])

    invisible(xy)
}


plotnd <- function(nMm, data=NULL,
                   diag.panel=NULL, 
               ...) {
    p <- nMm$dim

    # defining diagonal panel
    if (is.null(diag.panel)) {
        diag.panel <- function(i) {
            frame()
            box()
            text(0.5, 0.5, paste("var", i), cex=p)
        }
    }

    # setting up the plot
    opar <- par(mfcol=c(p,p) ) ## FIXME: some more: mar, oma
    on.exit(par(opar))


    # preallocating return data
    pts <- vector("list", p^2)
    dim(pts) <- c(p,p)

    for (i in 1:p) { 
        for (j in 1:p) {
            if (j == i) { # diagonal panel
                diag.panel(i)
            } else { # off diag panels, plotted with ellipsePts()
                pts[i,j][[1]] <- plot(reducednorMmix(nMm, c(i,j)), ...)
                if (!is.null(data)) points(data[,c(i,j)])
            }
        }
    }
    
    invisible(pts)
}


reducednorMmix <- function(nm, u) { 
    # nm: norMmix obj
    # u : dimension index
    #
    # returns: norMmix with reduced dimensions according to index vector u

    mu <- nm$mu[u,]
    Sig <- nm$Sigma[u,u,]
    return(norMmix(mu, Sigma=Sig, nm$weight, model=nm$model))
}


plot.norMmixMLE <- function(x, y=NULL, points=TRUE, ...) {
    plot(x$norMmix, data=if (points) x$x else NULL, ...)
    ## FIXME: sub-title or something about MLE (BIC, Likelihood, ..)
    # maybe something like title(sub="asdf"), might have to modify 
    # par()
}


# plot function for norMmix objects
# \code{plot.norMmix} returns invisibly coordinates of bounding ellipses of distribution
plot.norMmix <- function(x, y=NULL, ... ) {
    stopifnot(is.list(x), length(p <- x$dim) == 1)
    if (p == 2)
        plot2d(x, ...)
    else ## if (p>2)
        plotnd(x, ...)
}


### FIXME: Must document this ;  examples for *both* cases, ....
plot.fittednorMmix <- function(x, main="unnamed", plotbest=FALSE, ...) {
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

        matplot(bicmat, type="l", xlab="components", ylab="BIC", col=cl, lty=1:10, ...)
        legend("topright" , models, fill=cl, lty=1:10)
    } else {
        # TODO: like massplot (MM ??? )
        bk <- as.integer(best[1])
        bmodel <- best[2]
        plot(x$nMm[bk,bmodel][[1]]$norMmix, ...)
        points(x$x)
    }
    ## in both cases
    title(main=main)
    mtext(paste("best fit = ", best[1], best[2]))
}


