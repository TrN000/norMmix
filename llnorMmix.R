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
			      "VVI","EEE","VEE","EVV","VVV")
		      ) {

	# 1. sanity check on arguments
	# 2. transform par. to norMmix
	# 3. calculate log-lik
	# 4. return log-lik


	# 1. san check

	# 2. transform

	trafo <- match.arg(trafo)
	model <- match.arg(model)

	# 3. calc log-lik



	# get w

	w <- switch(trafo,
		    "clr1" = {if (k==1) 1
		    	     else clr1inv(par.[1:(k-1)])},

		    "logit" = logitinv(par.[1:(k-1)]),

		    stop("error in w switch in llnorMmix")
		    )

	if (!(sum(w)==1)) return(Inf)


	# get mu

	mu <- matrix(par.[k:(k+p*k-1)],p,k)


	# start of relevant parameters:

	f <- k + p*k # weights -1 + means +1 => start of alpha

	f1 <- f # end of alpha if uniform
	f2 <- f+k-1L # end of alpha if var

	f1.1 <- f1 +1L #start of D. if alpha unif.
	f2.1 <- f1 + k # start of D. if alpha varialbe

	f11 <- f1 + p # end of D. if D. uniform and alpha uniform
	f12 <- f1 + p*k # end D. if D. var and alpha unif.
	f21 <- f2 + p # end of D. if D. uniform and alpha variable
	f22 <- f2 + p*k # end of D. if D.var and alpha var

	f11.1 <- f11 +1L # start of L if alpha unif D unif
	f21.1 <- f21 +1L # start of L if alpha var D unif
	f12.1 <- f12 +1L # start of L if alpha unif D var
	f22.1 <- f22 +1L # start of L if alpha var D var

	f111 <- f11 + p*(p-1)/2 # end of L if alpha unif D unif
	f211 <- f21 + p*(p-1)/2 # end of L if alpha var D unif
	f121 <- f12 + k*p*(p-1)/2 # end of L if alpha unif D var
	f221 <- f22 + k*p*(p-1)/2 # end of L if alpha var D var


	retval <- 0

	retval <- switch(model,
		
	"EII" = {alpha <- par.[f]
		for (i in 1:k) {
			rss <- colSums(exp(alpha)*(t(x)-mu[,i])^2)
			retval <- retval+sum(-0.5*p*(alpha+log(2*pi))-0.5*rss)
		}
		retval},

	"VII" = {alpha <- par.[f:f2]
		for (i in 1:k) {
			rss <- colSums(exp(alpha[i])*(t(x)-mu[,i])^2)
			retval <- retval+sum(-0.5*p*(alpha[i]+log(2*pi))-0.5*rss)
		}
		retval},

	"EEI" = {alpha <- par.[f]
		D. <- par.[f1.1:f11]
		for (i in 1:k) {
			rss <- colSums(exp(alpha+D.)*(t(x)-mu[,i])^2)
			retval <- retval+sum(-0.5*p*(alpha+log(2*pi))-0.5*rss)
		}
		retval},

	"VEI" = {alpha <- par.[f:f2]
		D. <- par.[f2.1:f21]
		for (i in 1:k) {
			rss <- colSums(exp(alpha[i]+D.[,i])*(t(x)-mu[,i])^2)
			retval <- retval+sum(-0.5*p*(alpha[i]+log(2*pi))-0.5*rss)
		}
		retval},



	stop("error in temp switch in llnorMmix")
	)






#	y <- 0
#
#	for (i in 1:k) {
#		## this part only temporary until a faster solution is found
#		y <- y + w[i]*mvtnorm::dmvnorm(x,mean=mu[,i],sigma=sig[,,i])
#	}

	# 4. return


}
