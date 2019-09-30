### plotting methods for norMmix objects


plot2d.norMmix <- function(nMm, data, xlim=NULL, ylim=NULL, bounds=0.05,
                           type="l", lty=2, newWindow=TRUE, npoints=250,
                           col="red",  fill=TRUE, fillcolor="red",
	                   ... )
{
    w <- nMm$weight
    mu <- nMm$mu
    sig <- nMm$Sigma
    k <- nMm$k

    ## calculate smart values for xlim, ylim

    ellipsecoords <- vector()

    if ( is.null(xlim) || is.null(ylim) ) {
        for (i in 1:k) {
            ellipsecoords  <- rbind(ellipsecoords, mixtools::ellipse(mu=mu[,i], sigma=sig[,,i], newplot=FALSE, draw=FALSE, npoints=npoints))
        }
    }

    xbounds <- 0
    ybounds <- 0

    if (is.null(xlim)) {
        xlim <- c(min(ellipsecoords[,1]), max(ellipsecoords[,1]))
        xbounds <- bounds
    }

    if (is.null(ylim)) {
        ylim <- c(min(ellipsecoords[,2]), max(ellipsecoords[,2]))
        ybounds <- bounds
    }

    diffx <- abs(xlim[1]-xlim[2])*xbounds
    diffy <- abs(ylim[1]-ylim[2])*ybounds

    xlim <- c(xlim[1]-diffx, xlim[2]+diffx)
    ylim <- c(ylim[1]-diffy, ylim[2]+diffy)


    ## whether to plot new

    ifplot  <- vector("logical", k)

    if (newWindow) ifplot[1] <- TRUE

    ## determine fill color

    fco <- c(col2rgb(fillcolor)/255,(w[i]*0.8+0.1))

    ## add ellipses

    for (i in 1:k) {
        x <- mixtools::ellipse(mu=mu[,i], sigma=sig[,,i], newplot=ifplot[i], draw=TRUE, xlim=xlim, ylim=ylim,  type=type, lty=lty, col=col, npoints=npoints, ...)
        if (fill) polygon(x[,1], x[,2], col=rgb(red=fco[1],green=fco[2],blue=fco[3],alpha=fco[4]), border= NA )
    }

    ## label clusters

    text( mu[1,], mu[2,], sprintf("cluster %s", 1:k) )

    if (!is.null(data)) {
        points(data[,c(1,2)])
    }


    invisible(ellipsecoords)
}




plotnd.norMmix <- function(nMm,npoints=500, fillcolor="red",
                           alpha=0.05, ...)
{
    stopifnot( inherits(nMm, "norMmix") )

    w <- nMm$weight
    mu <- nMm$mu
    sig <- nMm$Sigma
    k <- nMm$k
    p <- nMm$dim


    npoints <- npoints*p


    ## get coords

    coord <- list()
    coarr <- matrix(0,k*npoints,p)
    corange <- list()

    for (i in 1:k) {
        r <- MASS::mvrnorm(n=npoints, mu=rep(0,p), sig[,,i])
        r <- apply(r,1, function(j) j/norm(j,"2"))
        r <- r*sqrt(qchisq(1-alpha,2))
        r <- sig[,,i]%*%r
        r <- r+mu[,i]
        r <- t(r)

        coord[[i]] <- r # cant use chull yet, only works on planar coords
        coarr[(1+(i-1)*npoints):(i*npoints),] <- r
        corange[[i]] <- apply(r,2,range)
    }


    ## color
    fco <- c(col2rgb(fillcolor)/255,(w[i]*0.8+0.1))
    fco <- rgb(red=fco[1],green=fco[2],blue=fco[3],alpha=fco[4])

    ploy <- function(x,y) {
        npoints <- eval.parent(npoints, n=2)
        fco <- eval.parent(fco, n=2)
        k <- eval.parent(k, n=2)
        w <- eval.parent(w, n=2)


        xs <- matrix(x,npoints,k)
        ys <- matrix(y,npoints,k)

        #points(x,y)
        for (i in 1:k) {
            ss <- cbind(xs[,i],ys[,i])
            polygon(ss[chull(ss),], col=fco)
        }
    }

    pairs(coarr, panel=ploy)

#    polygon(r[chull(r[,c(1,2)]),c(1,2)],col="red")


    invisible(coord)

}


#' plot function for norMmix objects
#'
#' \code{plot.norMmix} returns invisibly coordinates of bounding ellipses of distribution
#'
#' @export
plot.norMmix <- function(x, data=NULL, ... ) {
    ## TODO: make so data can also be missing
    stopifnot(is.list(x), length(p <- x$dim) == 1)
    if (p == 2)
        plot2d.norMmix(x, data, ... )
    else ## if (p>2)
        plotnd.norMmix(x, ...)
}





############################################################

plot.fittednorMmix <- function(x, name="unnamed", plotbest=FALSE, ...) {
    stopifnot(inherits(x, "fittednorMmix"))

    k <- x$k
    models <- x$models
    n <- x$n
    p <- x$p

    bicmat <- BIC(x)[[1]]
    best <- BIC(x)[[2]]

    cl <- rainbow(length(models))

    if (!plotbest) {
        matplot(bicmat, type="l", xlab="clusters", ylab="BIC", col=cl, lty=1:10, ...)
        title(main=name)
        legend("topright" , models, fill=cl, lty=1:10)
        mtext(paste("best fit = ", best[1], best[2]))
    } else {
        bk <- as.integer(best[1])
        bmodel <- best[2]
        plot(x$nMm[bk,bmodel][[1]]$norMmix, ...)
        title(main=name)
        mtext(paste("best fit = ", best[1], best[2]))
    }
}


