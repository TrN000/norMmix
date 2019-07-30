#### the extra m stands for multivariate




## norMmix constructor

tole <- 1000*.Machine$double.eps




hilf <- function(x, k){ # help function evals to true if x is sym pos def array
	for (i in 1:k){
		if ( !(isSymmetric(as.matrix(x[,,i]), tol=tole)&& matrixcalc::is.positive.semi.definite(x[,,i],tol=tole)) )
			stop("Sigma is not sym pos def")
	}

	return(TRUE)
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
	##	mu: matrix with columns as vector means of dist.
	##	Sigma: 	option 0: default, generates all identity
	##			covariance mat.
	##		option 1: scalar, generates EII dist
	##		option 2: vector of length k, generates VII
	##			distribution
	##		option 3: array of dimension p x p x k. 
	##			covariance matrices of distributions
	##	weight: vector of length k, sums to 1
	##	name: name attribute
	##	type: type of distribution VVV, IVV etc.
	## --------------------------------------------------------
	## Value: returns objext of class 'norMmix'
	## --------------------------------------------------------
	## Author: nicolas trutmann, Date:2019-06-19

	# p dimension 
	# k number of components

	# ispect mu 
	if (!is.numeric(mu)) stop("'mu' must be numeric")

	if ( is.vector(mu) ){
		k <- 1
		p <- length(mu)
		as.matrix(mu)
	} else if ( is.matrix(mu) ) {
		k <- ncol(mu)
		p <- nrow(mu)
	} else stop("mu is neither vector nor matrix")



	#ispect Sigma
	if (!is.numeric(Sigma)) stop("'Sigma' must be numeric")
	if (is.vector(Sigma) && length(Sigma)==1)
		Sigma <- array( diag(Sigma,p), c(p,p,k) )
	else if (is.vector(Sigma) && length(Sigma)==k)
		Sigma <- array(unlist(lapply(Sigma, function(j) diag(j,p))), c(p,p,k))
	else if (is.array(Sigma) && dim(Sigma) == c(p,p,k))
		Sigma <- Sigma
	else stop("'Sigma' not among recognized formats")
	if (!hilf(Sigma,k)) stop("error with sym pos def Sigma")

	#inspect weight
	if (!is.numeric(weight)) stop("'weight' must be numeric")
	if (! (is.vector(weight) && length(weight)==k) )
		stop("weight is not of correct dimension")
	if (!(all(weight >=0) && (sum(weight) - 1 < 1000*.Machine$double.eps)))
		stop("weight doesn't sum to 1 or isn't positive")

	if (missing(model)) {model <- "VVV"}
	else {model <- match.arg(model)}

	name <- sprintf("model = %s , clusters = %s", model, k)

	structure( name = name,
		  class = "norMmix",
		  .Data = list(mu = mu , Sigma = Sigma, weight = weight,
		  		k=k, dim=p, model=model)
	)

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

rnorMmix <- function(
		     n = 511,
		     obj
	  )
{

	## Purpose: generates random values distributed by NMD
	## -------------------------------------------------------------------
	## Arguments:
	## 	obj: of type norMmix
	##	n: number of values desired
       	## -------------------------------------------------------------------
	## Value: matrix, columns are vectors
	## -------------------------------------------------------------------
	## Author: nicolas trutmann, Date:2019-06-21

	if(!inherits(obj, "norMmix")) stop("'obj' must be of type norMmix")

	mu <- obj$mu
	Sigma <- obj$Sigma
	weight <- obj$weight
	p <- obj$dim

	#need case n=1 distiction here

	nj <- rmultinom(n=1, size=n, prob=weight)
#	a <- matrix(unlist(lapply( seq(along=nj), function(j)mvrnorm(n=nj[j], mu=mu[,j], Sigma=Sigma[,,j]) )), ncol=p, byrow=TRUE)
	## this approach doesnt work matrices arent concatenated properly
	a <- do.call( rbind,lapply( seq(along=nj), function(j) MASS::mvrnorm(n=nj[j], mu=mu[,j], Sigma=Sigma[,,j]) ))
}


dnorMmix <- function(x, nMm) {

	is.norMmix(nMm)

	p <- nMm$dim
	k <- nMm$k

	if (is.vector(x)) stopifnot(length(x)==p)
	if (is.matrix(x)) stopifnot(ncol(x)==p)

	ret <- 0

	for (i in 1:k) {
		ret <- ret + nMm$weight[i]*mvtnorm::dmvnorm(x, mean=nMm$mu[,i], sigma=nMm$Sigma[,,i])
	}

	ret
}



	






plot2d.norMmix <- function(nMm, xlim=NULL, ylim=NULL, bounds=0.05,
			   type="l", lty=2, newWindow=TRUE, npoints=250, 
			   col="red",  fill=TRUE, fillcolor="red")
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


	## if newWindow draw canvas from scratch
	
	if (newWindow) {
		plot.new()
		plot.window(xlim=xlim, ylim=ylim)
	}

	## determine fill color

	fco <- c(col2rgb(fillcolor)/255,(w[i]*0.8+0.1))

	## add ellipses 

	for (i in 1:k) {
		x <- mixtools::ellipse(mu=mu[,i], sigma=sig[,,i], newplot=FALSE, draw=TRUE, xlim=xlim, ylim=ylim,  type=type, lty=lty, col=col, npoints=npoints)
		if (fill) polygon(x[,1], x[,2], col=rgb(red=fco[1],green=fco[2],blue=fco[3],alpha=fco[4]), border= NA )
	}

	## axes and grid
		
	axis(1)
	axis(2)
	grid()

	## label clusters
	
	text( mu[1,], mu[2,], sprintf("cluster %s", 1:k) )


	invisible(ellipsecoords)
}


#' plot function for norMmix objects
#'
#' \code{plot.norMmix} returns invisibly coordinates of bounding ellipses of distribution
#'
#' This is the S3 method for plotting norMmix objects. atm only 2 dimensional objects are supported.
#'
#' @examples
#' plot(MW212) ; points(rnorMmix(n=500, MW212))
#' @export

plot.norMmix <- function(nMm, ... ) {

	if (nMm$dim==2) plot2d.norMmix(nMm, ... )

	else if (nMm$dim>2) warning("methods for more than 2 dim not yet done")


}



metric.norMmix <- function(n1,n2, type="2", matchby=c("mu","id")) {

	stopifnot( is.norMmix(n1), is.norMmix(n2) )
	stopifnot( isTRUE(all.equal(n1$k, n2$k)) )

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
