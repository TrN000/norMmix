#### the extra m stands for multivariate

## norMmix constructor





## Auxiliary function evals to TRUE if x is sym pos def array, otherwise character string with msg
okSigma <- function(Sig, tol1 = 1000*.Machine$double.eps, tol2 = 1e-10) {
    if(!is.numeric(Sig) || !is.array(Sig) || length(d <- dim(Sig)) != 3)
        "is not a numeric array of rank 3"
    else if(prod(d) == 0)
        "at least one of the dimensions is zero"
    else { ## all d[.] >= 1
        k <- d[[3]]
        p <- d[[1]]
        if(d[[2]] != p)
            "Sigma matrix dimension is not square"
        else {
            for (i in 1:k) {
                Si <- Sig[,,i]
                if(!isSymmetric(Si, tol=tol1))
                    return(paste0("Sigma[ , ,",i,"] is not symmetric"))
                ## else
                ev <- eigen(Si, only.values=TRUE)$values
                if(any(ev < -tol2))
                    return(paste0("Sigma[ , ,",i,"] is not positive semi-definite"))
            }
            TRUE
        }
    }
}





#' Constructor for nMm 'objects'
#'
#' \code{norMmix} returns structure containing defining parameters of a normal
#' mixture
#'
#'
#'
#' @param mu matrix of means. should mu be a vector it will assume k=1
#' to circumvent this behavoiur use as.matrix(mu) beforehand
#' @param Sigma array of covariance matrices
#' @param weight weights of mixture model components
#' @param name gives the option of naming mixture
#' @param model see desc
#'
#' @export

norMmix <- function(
            mu,
            Sigma = NULL,
            weight = rep(1/k, k),
            name = NULL,
            model= c("EII","VII","EEI","VEI","EVI",
                     "VVI","EEE","VEE","EVV","VVV")
            )
{
    ## Purpose: constructor for 'norMmix' (multivariate normix)
    ## makes sure values are sane.
    ## --------------------------------------------------------
    ## Arguments:
    ##    mu: matrix with columns as vector means of dist.
    ##    Sigma:     option 0: default, generates all identity
    ##            covariance mat.
    ##        option 1: scalar, generates EII dist
    ##        option 2: vector of length k, generates VII
    ##            distribution
    ##        option 3: array of dimension p x p x k.
    ##            covariance matrices of distributions
    ##    weight: vector of length k, sums to 1
    ##    name: name attribute
    ##    type: type of distribution VVV, IVV etc.
    ## --------------------------------------------------------
    ## Value: returns objext of class 'norMmix'
    ## --------------------------------------------------------
    ## Author: nicolas trutmann, Date:2019-06-19

    stopifnot(is.numeric(mu))
    if(!is.matrix(mu)) mu <- as.matrix(mu) # p x 1  typically
    p <- nrow(mu) # p = dimension
    k <- ncol(mu) # k = number of components

    dS <- c(p,p,k) # == dim(Sigma)
    if(is.null(Sigma)) Sigma <- 1  # (option 0)
    else stopifnot(is.numeric(Sigma))
    isArr <- is.array(Sigma) # if not array, is also not matrix
    if (!isArr && length(Sigma) == 1)
        Sigma <- array(diag(Sigma,p), dS)
    else if(!isArr && length(Sigma) == k)
        Sigma <- array(vapply(Sigma, function(s) diag(s,p), diag(p)), dS)
    else if(!(isArr && all(dim(Sigma) == c(p,p,k))))
        stop("'Sigma' not among recognized formats")
    if (!isTRUE(m <- okSigma(Sigma))) stop(m)

    # inspect weight
    stopifnot(is.numeric(weight))
    if(length(weight) != k)
        stop("weight is not of length k")
    if (!(all(weight >= 0) && (abs(sum(weight) - 1) < 1000*.Machine$double.eps)))
        stop("weight doesn't sum to 1 or isn't positive")

    model <- if(missing(model)) "VVV" else match.arg(model)
    if(is.null(name))
        name <- sprintf("model \"%s\", G = %s", model, k)
    structure(name = name,
              class = "norMmix",
              list(model = model,
                   mu = mu, Sigma = Sigma, weight = weight
                 , k = k  # == length(weight) == ncol(mu)
                 , dim = p # == nrow(mu)
                   ) )
}






is.norMmix <- function(obj){
    inherits(obj, "norMmix")
}






mean.norMmix <- function(obj){
    if (!is.norMmix(obj)) stop("object is not norMmix")
    k <- obj$k
    mu <- obj$mu
    w <- obj$weight

    me <- rep(0, obj$dim)

    for (i in 1:k) {
        me <- me + w[i]*mu[,i]
    }
    return(me)
}






### rnorMmix

rnorMmix <- function(n, obj, index=FALSE, permute=TRUE)
{

    ## Purpose: generates random values distributed by NMD
    ## -------------------------------------------------------------------
    ## Arguments:
    ##    n: number of p-dimensional observations desired
    ##  obj: of type norMmix
    ## -------------------------------------------------------------------
    ## Value: matrix  n x p (columns are vectors)
    ## -------------------------------------------------------------------
    ## Author: nicolas trutmann, Date:2019-06-21

    if(!inherits(obj, "norMmix")) stop("argument must be of class \"norMmix\"")

    mu <- obj$mu
    Sigma <- obj$Sigma
    weight <- obj$weight
    p <- obj$dim
    nj <- rmultinom(n=1, size=n, prob = weight)
### FIXME: again this does  chol(Sigma_j)  j = 1...k .. when actually we could've rather stored  chol(Sigma) instead of Sigma
    a <- do.call(rbind,
                 lapply(seq_along(nj), function(j) mvrnorm(n=nj[j], mu=mu[,j], Sigma=Sigma[,,j])))
    if (index) {
        cl <- rep(seq_along(weight), times=nj)
        a <- cbind(cl, a)
    }
    if(permute)
        a <- a[sample.int(n),]
    else
        a
}


dnorMmix <- function(x, nMm) {
    stopifnot(is.norMmix(nMm), is.numeric(x),
              length(p <- nMm$dim) == 1, p >= 1,
              length(k <- nMm$ k ) == 1, k >= 1)
    if(!is.matrix(x)) x <- as.matrix(x)
    stopifnot(ncol(x) == p)
    ## FIXME: Using dmvnorm() is slow as it needs to solve( Sigma_j ) for each  j = 1..k
    ret <- 0
    for (i in 1:k) {
        ret <- ret + nMm$weight[i]*mvtnorm::dmvnorm(x, mean=nMm$mu[,i], sigma=nMm$Sigma[,,i])
    }
    ret
}








plot2d.norMmix <- function(nMm, xlim=NULL, ylim=NULL, bounds=0.05,
               type="l", lty=2, newWindow=TRUE, npoints=250,
               col="red",  fill=TRUE, fillcolor="red",
	           ...) {
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


    invisible(ellipsecoords)
}

plotnd.norMmix <- function(nMm,npoints=500, fillcolor="red",
                           alpha=0.05, ...) {
    stopifnot( inherits(nMm, "norMmix") )

    w <- nMm$weight
    mu <- nMm$mu
    sig <- nMm$Sigma
    k <- nMm$k
    p <- nMm$dim


    npoints <- npoints*p


    ## get coords??

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
plot.norMmix <- function(x, ... ) {
    stopifnot(is.list(x), length(p <- x$dim) == 1)
    if (p == 2)
        plot2d.norMmix(x, ... )
    else ## if (p>2)
        plotnd.norMmix(x, ...)
}



metric.norMmix <- function(n1,n2, type="2", matchby=c("mu","id")) {

    stopifnot( is.norMmix(n1), is.norMmix(n2) )
    stopifnot( all.equal(n1$k, n2$k) )

    matchby <- match.arg(matchby)

    k <- n1$k

    # sort cluster to compare by difference in means

    order. <- switch(matchby,

        "id" = 1:k,

        "mu" = {
                order. <- vector()
                m1 <- n1$mu
                m2 <- n2$mu
                for (i in 1:k) {
                    diffmu <- apply((m1-m2[,i]),2,norm,type=type)
                    order.[i] <-  which.min(diffmu)
                }
                order.
            },

        stop("error in matchby statement")
        )


    deltamu <- apply(m1-m2[,order.],2,norm,type=type)

    deltasig <- apply(n1$Sigma-n2$Sigma[,,order.],3,norm,type=type)

    deltaweight <- n1$weight - n2$weight[order.]

    ## some penalty value

    p <- n1$dim

    pmu <- sqrt(deltamu^2)/(p*k)
    psig <- sqrt(deltasig^2)/(p*p*k)
    pweight <- sqrt(deltaweight^2)/k

    penalty <- sum( pmu+psig+pweight )

    list(order.=order., deltamu=deltamu, deltasig=deltasig, deltaweight=deltaweight, penalty=penalty)

}
