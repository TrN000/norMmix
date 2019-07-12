#### MLE for norMmix objects

#' @include norMmixEM.R

#' Maximum likelihood Estimation for normal mixture models
#' 
#' \code{norMmixMLE} returns fitted nMm obj
#' 
#' @param x
#' @param trafo
#' @param model
#' 
#' @export

norMmixMLE <- function(
		       x,
		       k,
		       trafo = c("clr1", "logit"),
		       model = c("EII","VII","EEI","VEI",
			      "EVI","VVI","EVV","VVV")
		       ) {
	
	# 1. san check call
	# 2. prep nMm obj
	# 3. apply optim
	# 4. return


	# 1.

	trafo <- match.arg(trafo)
	model <- match.arg(model)

	stopifnot(is.numeric(x), is.numeric(k), (n <- ncol(x))>1)
	p <- nrow(x)

	k <- as.integer(k)



	# 2.

	par.temp <- mstep.nMm(x, tau, mu, Sigma, weight)


	# 3.

	neglogl <- function(par.) -llnorMmix(par., x=x, p=p, k=k, trafo=trafo, model=model)


	#optr <- #here optidd


	# 4.

	list(weight=weight, mu=mu, Sigma=Sigma)

}
