#### the extra m stands for multivariate




## norMmix constructor

tole <- 1000*.Machine$double.eps

# do I really need this many?
library(abind) # can replace these with base function

Shilf <- function(L, p){ # help function for norMmix
	Shilfarray <- array(0, c(p,p,length(L)))
	for (i in 1:length(L)){
		Shilfarray[,,i] <- L[i]*diag(p)
	}
	return(Shilfarray)
}








hilf <- function(x, k){ # help function evals to true if x is sym pos def array
	for (i in 1:k){
		if ( !(isSymmetric(x[,,i], tol=tole)&& matrixcalc::is.positive.semi.definite(x[,,i],tol=tole)) )
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
#' @param mu matrix of means
#' @param Sigma array of covariance matrices
#' @param weight weights of mixture model components
#' @param name gives the option of naming mixture
#' @param model see desc
#'
#' @export

norMmix <- function(
		    mu,
		    Sigma = abind( rep(list(diag(p)), k), along = 3),
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
	if (is.matrix(mu) == FALSE){
		p <- 1
		k <- length(mu)
	} else {
		k <- ncol(mu)
		p <- nrow(mu)
	}


	#ispect Sigma
	if (!is.numeric(Sigma)) stop("'Sigma' must be numeric")
	if (is.vector(Sigma) && length(Sigma)==1)
		Sigma <- abind( rep(list(Sigma*diag(p)), k), along = 3)
	else if (is.vector(Sigma) && length(Sigma)==k)
		Sigma <- Shilf(L=Sigma, p=p) 
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

	structure( name = name,
		  class = "norMmix",
		  .Data = list(mu = mu , Sigma = Sigma, weight = weight,
		  		k=k, dim=p, model=model)
	)

	#nMm <- list(mu = mu , Sigma = Sigma, weight = weight,
	#	    k=k, dim=p, model=model)
        #
	#attr(nMm,"name") <- name
        #
	#class(nMm) <- "norMmix"



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
		     obj,
		     n = 511
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
	a <- do.call( rbind,lapply( seq(along=nj),function(j)MASS::mvrnorm(n=nj[j], mu=mu[,j], Sigma=Sigma[,,j]) ))
}
