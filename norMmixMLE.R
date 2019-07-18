#### MLE for norMmix objects

#' @include norMmixEM.R





{}


#' Maximum likelihood Estimation for normal mixture models
#' 
#' \code{norMmixMLE} returns fitted nMm obj
#' 
#' Uses clara() and one M-step from EM-algorithm to initialize parameters
#' after that uses general optimizer optim() to calculate ML.
#'
#' @param x sample matrix
#' @param k number of clusters
#' @param trafo transformation to be used
#' @param model model to be assumed
#' 
#' @export
norMmixMLE <- function(
		       x,
		       k,
		       trafo = c("clr1", "logit"),
		       model = c("EII","VII","EEI","VEI","EVI",
				 "VVI","EEE","VEE","EVV","VVV")
		       ) {
	
	# 1. san check call
	# 2. prep nMm obj
	# 3. apply optim
	# 4. return


	# 1.

	trafo <- match.arg(trafo)
	model <- match.arg(model)

	stopifnot(is.numeric(x), is.numeric(k), (n <- nrow(x))>1)
	p <- ncol(x)

	k <- as.integer(k)

	#init tau??

	clus <- cluster::clara(x=x, k)

	index <- clus$clustering

	tau <- matrix(0,n,k)
	tau[cbind(1:n,index)] <- 1


	# 2.

	# one m-step
	nMm.temp <- mstep.nMm(x, tau, mu, Sigma, weight, k, p)

	# create par. vector out of m-step
	par. <- nMm2par(obj=nMm.temp, trafo=trafo, model=model)


	# 3.

	# define function to optimize as negative log-lik
	# also reduces the number of arguments to par.
	neglogl <- function(par.) {
		-llnorMmix(par.,x=x,p=p,k=k,trafo=trafo,model=model)
		}

	optr <- optim(par., neglogl, method = "BFGS")


	# 4.

	nMm <- par2nMmMLE(optr$par, p, k, trafo=trafo, model=model)

}
