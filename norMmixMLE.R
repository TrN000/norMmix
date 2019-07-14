#### MLE for norMmix objects

#' @include norMmixEM.R





{}


#' Maximum likelihood Estimation for normal mixture models
#' 
#' \code{norMmixMLE} returns fitted nMm obj
#' 
#'
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
	par.temp <- nMm2par(obj=nMm.temp, trafo=trafo, model=model)

	# log of alpha, D. 
	par. <- par2nMmMLE_inv(par.temp, p, k, trafo=trafo, model=model)


	# 3.

	neglogl <- function(par.) {
		-llnorMmix(par.,x=x,p=p,k=k,trafo=trafo,model=model)
		}

	optr <- optim(par., neglogl, method = "BFGS")

	nMm <- par2nMmMLE(optr$par, p, k, trafo=trafo, model=model)


	# 4.

	nMm
}
