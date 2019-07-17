#### the llnorMmix function, calculating log likelihood for a given
#### parameter vector

## Author: Nicolas Trutmann 2019-07-06

#' @include param.R




{}



#' Log-likelihood of parameter vector given data
#'
#' \code{llnorMmix} returns scalar log-likelihood 
#'
#'
#' description
#'
#' @param par. parameter vector
#' @param x sample matrix
#' @param p dimension
#' @param k number of clusters
#' @param trafo either centered log ratio or logit
#' @param model assumed distribution model of normal mixture
#'
#' @export

llnorMmix <- function(par., x, p, k,
		      trafo=c("clr1", "logit"),
		      model=c("EII","VII","EEI","VEI","EVI",
			      "VVI","EEE","VEE","EVV","VVV"),
		      MLE=FALSE ) {

	# 1. sanity check on arguments
	# 2. transform par. to norMmix
	# 3. calculate log-lik
	# 4. return log-lik


	# 1. san check

	# 2. transform

	trafo <- match.arg(trafo)
	model <- match.arg(model)

	if (MLE) {nMmobj <- par2nMm(par., p, k, trafo=trafo, model=model)}
	else {nMmobj <- par2nMmMLE(par., p, k, trafo=trafo, model=model)}

	# 3. calc log-lik

	w <- nMmobj$w
	mu <- nMmobj$mu
	sig <- nMmobj$Sig

	y <- 0

	for (i in 1:k) {
		y <- y + w[i]*mvtnorm::dmvnorm(x,mean=mu[,i],sigma=sig[,,i])
	}

	# 4. return

	res <- sum(log(y))

}
