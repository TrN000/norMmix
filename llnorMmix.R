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

	stopifnot( (ncol(x)==p || nrow(x)==p) )
	if (ncol(x)==p && nrow(x)!=p) x <- t(x)

	# 2. transform

	trafo <- match.arg(trafo)
	model <- match.arg(model)


	l2pi <- log(2*pi)

	# 3. calc log-lik



	# get w

	w <- switch(trafo,
		    "clr1" = {if (k==1) 1
		    	     else clr1inv(par.[1:(k-1)])},

		    "logit" = logitinv(par.[1:(k-1)]),

		    stop("error in w switch in llnorMmix")
		    )

	if (!(sum(w)==1)) return(-Inf)


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
			rss <- colSums((1/exp(alpha))*(x-mu[,i])^2)
			# this is vector of length n=sample size
			# calculates (x-mu)t * Sigma^-1 * (x-mu) for diagonal
			# cases.
			retval <- retval+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
			# adds likelihood of one component to retval
			# the formula in exp() is the log of likelihood
			# still of length n,
			# here we implicitly sum over components
		}
		sum(log(retval))},
		# here transform back to log
		# then sum over sample

	"VII" = {alpha <- par.[f:f2]
		for (i in 1:k) {
			rss <- colSums((1/exp(alpha[i]))*(x-mu[,i])^2)
			retval <- retval+w[i]*exp(-0.5*p*(alpha[i]+l2pi)-0.5*rss)
		}
		sum(log(retval))},

	"EEI" = {alpha <- par.[f]
		D. <- par.[f1.1:f11]
		for (i in 1:k) {
			rss <- colSums((1/exp(alpha+D.))*(x-mu[,i])^2)
			retval <- retval+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
		}
		sum(log(retval))},

	"VEI" = {alpha <- par.[f:f2]
		D. <- par.[f2.1:f21]
		for (i in 1:k) {
			rss <- colSums((1/exp(alpha[i]+D.))*(x-mu[,i])^2)
			retval <- retval+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
		}
		sum(log(retval))},

	"EVI" = {alpha <- par.[f]
		D. <- matrix(par.[f1.1:f12],p,k)
		for (i in 1:k) {
			rss <- colSums((1/exp(alpha+D.[,i]))*(x-mu[,i])^2)
			retval <- retval+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
		}
		sum(log(retval))},

	"VVI" = {alpha <- par.[f:f2]
		D. <- matrix(par.[f2.1:f22],p,k)
		for (i in 1:k) {
			rss <- colSums((1/exp(alpha[i]+D.[,i]))*(x-mu[,i])^2)
			retval <- retval+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
		}
		sum(log(retval))},

	"EEE" = {alpha <- par.[f]
		D. <- par.[f1.1:f11]
		L. <- diag(1,p)
		L.[lower.tri(L., diag=FALSE)] <- par.[f11.1:f111]
		for (i in 1:k) {
			rss <- colSums((1/exp(alpha+D.))*backsolve(L.,(x-mu[,i]))^2)
			retval <- retval+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
			print(sum(log(retval)))
		}
		sum(log(retval))},

	"VEE" = {alpha <- par.[f:f2]
		D. <- par.[f2.1:f21]
		L. <- diag(1,p)
		L.[lower.tri(L., diag=FALSE)] <- par.[f21.1:f211]
		for (i in 1:k) {
			rss <- colSums((1/exp(alpha[i]+D.))*backsolve(L., (x-mu[,i]))^2)
			retval <- retval+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
		}
		sum(log(retval))},

	"EVV" = {alpha <- par.[f]
		D. <- matrix(par.[f1.1:f12],p,k)
		L.temp <- matrix(par.[f12.1:f121],p*(p-1),k)
		for (i in 1:k) {
			L. <- diag(1,p)
			L.[lower.tri(L., diag=FALSE)] <- L.temp[,i]
			rss <- colSums((1/exp(alpha+D.[i]))*backsolve(L., (x-mu[,i]))^2)
			retval <- retval+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
		}
		sum(log(retval))},

	"VVV" = {alpha <- par.[f:f2]
		D. <- matrix(par.[f2.1:f22],p,k)
		L.temp <- matrix(par.[f22.1:f221],p*(p-1),k)
		for (i in 1:k) {
			L. <- diag(1,p)
			L.[lower.tri(L., diag=FALSE)] <- L.temp[,i]
			rss <- colSums((1/exp(alpha[i]+D.[,i]))*backsolve(L., (x-mu[,i]))^2)
			retval <- retval+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
		}
		sum(log(retval))},



	stop("error in temp switch in llnorMmix")
	)


retval




}





llmvtnorm <- function(par., x, p, k,
		      trafo=c("clr1","logit"),
		      model=c("EII","VII","EEI","VEI","EEE",
			      "VEE","EVI","VVI","EVV","VVV")
		      ) {


	p <- as.integer(p)
	k <- as.integer(k)

	tr <- match.arg(trafo)
	mo <- match.arg(model)

	if (ncol(x)!=p && nrow(x)==p) x <- t(x)



	nmm <- par2nMm(par., p, k, trafo=tr, model=mo)

	w <- nmm$w
	mu <- nmm$mu
	sig <- nmm$Sigma


	y <- 0

	for (i in 1:k) {
		y <- y + w[i]*mvtnorm::dmvnorm(x,mean=mu[,i],sigma=sig[,,i])
		print(sum(log(y)))
	}


	sum(log(y))

}
