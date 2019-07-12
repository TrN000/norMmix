#### the llnorMmix function, calculating log likelihood for a given
#### parameter vector

## Author: Nicolas Trutmann 2019-07-06

#' @include param.R



#' log likelihood of parameter vector given data
#' 

llnorMmix <- function(par., x, p, k,
		      trafo=c("clr1", "logit"),
		      model=c("EII","VII","EEI","VEI",
			      "EVI","VVI","EVV","VVV") ){

	## Purpose: calculates log-likelihood of fitted norMmix model
	## ------------------------------------------------------------
	## Arguments: par.: parameter vector
	##	p: dimension
	##	k: number of clusters
	##	trafo: either centered log ratio or logit
	##	model: assumed distribution model of normal mixture
	## ------------------------------------------------------------
	## value: integer
	## ------------------------------------------------------------
	## Author: Nicolas Trutmann
	# 1. sanity check on arguments
	# 2. transform par. to norMmix
	# 3. calculate log-lik
	# 4. return log-lik


	# 1. san check

	# 2. transform

	nMmobj <- par2nMm(par., p, k,
			  trafo=c("clr1", "logit"),
			  model=c("EII","VII","EEI","VEI",
				  "EVI","VVI","EVV","VVV"),
			  MLE=TRUE) 

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
	res

}
