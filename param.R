#### functions handling parameter manipulation. par2nM nM2par etc

### for info on models see Celeux and Govaert 1995

#' @include Cholesky.R
#' @include weight.R


## map lower.tri to vec
ld. <- function(mat){
	x <- mat[lower.tri(mat, diag=FALSE)]
	x
}

## map vec to lower.tri
dl. <- function(d,x,p){
	mat <- matrix(0,p,p)
	mat[lower.tri(mat,diag=FALSE)==TRUE] <- x
	mat.temp <- mat + diag(rep(1,p))  
	mat <- mat.temp %*% diag(d) %*% t(mat.temp)
	mat
}


#' normal multivariate mixture model to parameter for MLE
#'
#' \code{nMm2par} returns vector of parameters of norMmix objects
#'
#' This transformation forms a vector from the parameters of a normal
#' mixture. These consist of weights, means and covariance matrices.
#' Weights are transformed according to 'trafo' param; means are 
#' unchanged.
#' Cov mats are given as D and L from the LDLt decomposition
#'
#' @param obj list containing sig= covariance matrix array, mu= mean vector matrix, w= weights, k= number of clusters, p= dimension
#' @param trafo either "clr1" or "logit"
#' @param model one of "asdf..."
#' @examples
#' A  <- MW2nm4
#' nMm2par( A, trafo="clr1", model=A$model )
#'
#' @export

nMm2par <- function(obj, 
		    trafo=c("clr1", "logit"),
		    model=c("EII","VII","EEI","VEI","EVI",
			    "VVI","EEE","VEE","EVV","VVV")
		    ){

	#transferring values of obj to handier variables
	w <- obj$w
	mu <- obj$mu
	sig <- obj$Sigma
	p <- obj$dim
	k <- obj$k

	trafo <- match.arg(trafo)
	model <- match.arg(model)

	##checks

	# weights

	stopifnot( isTRUE(all.equal(sum(w),1)), (length(w)==k), is.numeric(w) )

	if ( any(is.na(w))||any(is.nan(w))||
	     any(is.null(w))||any(is.infinite(w)) ) {
		warning("result may contain: NA, Nan, NULL, Inf in weight")
	}

	# mu

	stopifnot( (dim(mu)==c(p,k)), (is.numeric(mu)) )

	if ( any(is.na(mu))||any(is.nan(mu))||
	     any(is.null(mu))||any(is.infinite(mu)) ) {
		warning("result may contain: NA, Nan, NULL, Inf in mu")
	}

	# Sigma

	stopifnot( (dim(sig)==c(p,p,k)), (is.numeric(sig)) )
	stopifnot( isTRUE( all( apply(sig,3, function(j) (ldl(j)$Diag >= 0 )))))

	if ( any(is.na(sig))||any(is.nan(sig))||
	     any(is.null(sig))||any(is.infinite(sig)) ) {
		warning("result may contain: NA, Nan, NULL, Inf in sig")
	}


	#output vector of parameter values
	c(
	  w <- switch(trafo, #weights either logit or centered log ratio
		"logit" = logit(w),

    		"clr1" = clr1(w),

		stop("Error in weight trafo, ",trafo)
	    	),
	  mu, #means
	  Sigma <- switch(model, #model dependent covariance values
		"EII" = sig[1,1,1],

		"VII" = sig[1,1,],

		"EEI" = {D.temp <- diag(sig[,,1])
			alpha <- (prod(D.temp)^(1/p))
			D. <- D.temp/alpha
			c(alpha, D.)
			},

		"VEI" = {alpha <- apply(sig, 3, function(j) prod(diag(j))^(1/p) )
			D. <- diag(sig[,,1])/alpha[1]
			c(alpha, D.)
			},
			  
		"EVI" = {alpha <- prod(diag(sig[,,1]))^(1/p)
			D. <- apply(sig,3,diag)/alpha
			c(alpha,D.)
			},

		"VVI" = {alpha <- apply(sig,3, function(j) det(j)^(1/p)) 
			D.temp <- apply(sig,3,diag)
			D. <- D.temp %*% diag(1/alpha) # this is fastest?? https://stackoverflow.com/questions/20596433/how-to-divide-each-row-of-a-matrix-by-elements-of-a-vector-in-r
			c(alpha,D.)
			},

		"EEE" = {alpha <- prod( ldl(sig[,,1])$Diag )^(1/p)
			A. <- ldl(sig[,,1])
			c(alpha, A.$Diag/alpha, ld.(A.$L))},

		"VEE" = {alpha <- apply(sig,3, function(j) prod(ldl(j)$Diag )^(1/p))
			A. <- ldl(sig[,,1])
			c(alpha, A.$Diag/alpha[1], ld.(A.$L))},

		"EVV" = {alpha <- prod( ldl(sig[,,1])$Diag )^(1/p)
			D. <- apply(sig,3, function(j) ldl(j)$Diag/alpha)
			L. <- apply(sig,3, function(j) ld.( ldl(j)$L ))
			c(alpha, D., L.)},

		"VVV" = {alpha <- apply(sig,3, function(j) prod(ldl(j)$Diag )^(1/p))
			D.temp <- apply(sig,3, function(j) ldl(j)$Diag)
			D. <- D.temp %*% diag(1/alpha)
			L. <- apply(sig,3, function(j) ld.( ldl(j)$L ))
			c(alpha, D., L.)},

		stop("invalid argument in 'model'")
		)
	)
}





#' transform of parameter vector to normal mixture
#'
#' \code{par2nMm} returns list containing weight, mu, Sigma, k, dim
#'
#' this is the inverse function to nMm2par. Given a numeric vector
#' dimension and cluster number this function reconstructs a normal mixture
#' object.
#'
#' @param par. numeric vector of parameters
#' @param p dimension of space
#' @param trafo either "clr1" or "logit"
#' @param model See description 
#' 
#' @return returns this list: list(weight=w, mu=mu, Sigma=Sigma, k=k, dim=p)
#' @export

par2nMm <- function(par., p, k,
		    trafo=c("clr1", "logit"),
		    model=c("EII","VII","EEI","VEI","EEE",
			    "VEE","EVI","VVI","EVV","VVV")
		    ) {

	trafo <- match.arg(trafo)
	model <- match.arg(model)

	f <- k+p*k # start of Sigma/alpha
	f1 <- f + k - 1L # end of alpha if variable

	if (k ==1) {w.temp <- vector()}
	else {w.temp <- par.[1:(k-1)]}
	w <- switch(trafo,
	"logit" = {w.temp <- plogis(w.temp)
		   if ((sp. <- sum(w.temp)) > 1) {
			stop("weights add to greater than 1")
		   }
		 
     		   w <- c((1-sp.), w.temp)	},

	"clr1" = clr1inv(w.temp), 
	)

	mu <- matrix(par.[k:(k+p*k-1)], p, k)



	Sigma <- switch(model,

	# diagonal cases
	"EII" = {lambda <- par.[f]

		 sig <- array( rep(diag(lambda, p),k), c(p,p,k) )},

	"VII" = {lambda <- par.[f:f1]

		 sig <- array(unlist(lapply( lambda, function(j) diag(j,p) )), c(p,p,k)) },

	"EEI" = {lambda <- par.[f]

		 D.temp <- par.[f+1:f+p]

		 sig <- array( rep(diag(lambda*D.temp),k), c(p,p,k) )},

	"VEI" = {lambda <- par.[f:f1]

		 D.temp <- par.[f1+1:f1+p]

		 D. <- tcrossprod(D.temp,lambda)

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	"EVI" = {lambda <- par.[f]

		 D.temp <- matrix(par.[f1+1:f1+p*k],p,k)

		 D. <- D.temp*lambda

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	"VVI" = {lambda <- par.[f:f1]

		 D.temp <- matrix(par.[f1+1:f1+p*k],p,k)

		 D. <- D.temp %*% diag(lambda)

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	# variable cases

	"EEE" = {lambda <- par.[f]
		 D.temp <- par.[f+1:f+p]
		 D. <- D.temp*lambda

		 f3 <- (p*(p-1)/2)

		 L. <- par.[f+p+1:f+p+f3]
		 A. <- dl.(D.,L.,p)

		 sig.temp <- array(0, c(p,p,k))
		 for (i in i:k){
			 sig[,,i] <- A.
		 }
		 sig},


	"VEE" = {lambda <- par.[f:f1]
		 D.temp <- par.[f1+1:f1+p]

		 D. <- tcrossprod(D.temp,lambda)

		 f3 <- (p*(p-1)/2)

		 L. <- par.[f1+p+1:f1+p+f3]

		 sig.temp <- array(0, c(p,p,k))
		 for (i in i:k){
			 sig[,,i] <- dl.(D.[,i],L.,p)
		 }
		 sig},

	"EVV" = {lambda <- par.[f]

		 D.temp <- matrix(par.[f1+1:f1+p*k],p,k)

		 f2 <- f1 + p*k +1 # start of L_k
		 f3 <- (p*(p-1)/2)
		 f4 <- (p*(p-1)/2)*k
		 f5 <- f2 + f4 -1

		 L.temp <- matrix(par.[f2:f5],f3,k)

		 sig <- array(0, c(p,p,k))
		 for (i in 1:k) {
			 sig[,,i] <- dl.(D.[,i],L.temp[,i],p)
		 }
		 sig},

	"VVV" = {lambda <- par.[f:f1]

		 D.temp <- matrix(par.[(f1+1):(f1+p*k)],p,k)

		 D. <- D.temp %*% diag(lambda)

		 f2 <- f1 + p*k +1 # start of L_k
		 f3 <- (p*(p-1)/2)
		 f4 <- (p*(p-1)/2)*k
		 f5 <- f2 + f4 -1

		 L.temp <- matrix(par.[f2:f5],f3,k)

		 sig <- array(0, c(p,p,k))
		 for (i in 1:k) {
			 sig[,,i] <- dl.(D.[,i],L.temp[,i],p)
		 }
		 sig},
	stop("error in Sigma switch statement")
	)


	list( weight=w, mu=mu, Sigma=Sigma, k=k, dim=p )

}




#' wrapper function for par2nMm that preprocesses for norMmixMLE
#' 
#' \code{par2nMmMLE} returns nMm list guarantees that Sigma are sym. pos def
#' and alpha are nonnegative.
#'
#' This variant of par2nMm applies exp to mixture parameters alpha and D.
#' to guarantee that Sigma are always sym pos def. Note here that this does not
#' affect weight parameter. Those may still be out of bounds of permissible 
#' values.
#' 
#' @inheritParams par2nMm
#'
#' @export

par2nMmMLE <- function(par., p, k,
		    trafo=c("clr1", "logit"),
		    model=c("EII","VII","EEI","VEI","EEE",
			    "VEE","EVI","VVI","EVV","VVV")
		    ) {

	# start of relevant parameters:

	f <- k + p*k # weights -1 + means +1 => start of alpha
	f1 <- f # end of alpha if alpha uniform
	f1.1 <- f1 +1L #start of D. if alpha unif.
	f1.2 <- f1 + p # end of D. if D. uniform and alpha uniform
	f1.3 <- f1 + p*k # end D. if D. var and alpha unif.
	
	f2.1 <- f1 + k # start of D. if alpha varialbe
	f2.2 <- f2.1 + p # end of D. if D. uniform and alpha variable
	f2.3 <- f2.1 + p*k # end of D. if D.var and alpha var

	#only important ones are f1.2, f1.3, f2.2, f2.3


	model <- match.arg(model)

	par. <- switch(model,

	"EII" = {par.[f] <- exp(par.[f])
		 par.},

	"VII" = {par.[f:f2] <- exp(par.[f:f2])
		 par.},

	"EEI" = ,
	"EEE" = {par.[f:f1.1] <- exp(par.[f:f1.1])
		 par.},

	"VEI" = ,
	"VEE" = {par.[f:f2.2] <- exp(par.[f:f2.2])
		 par.},

	"EVI" = ,
	"EVV" = {par.[f:f1.3] <- exp(par.[f:f1.3])
		 par.},

	"VVI" = ,
	"VVV" = {par.[f:f2.3] <- exp(par.[f:f2.3])
		 par.}
	)


	par2nMm(par., p, k, trafo=trafo, model=model)

}






par2nMmMLE_inv <- function(par., p, k,
		    trafo=c("clr1", "logit"),
		    model=c("EII","VII","EEI","VEI","EEE",
			    "VEE","EVI","VVI","EVV","VVV")
		    ) {

	# start of relevant parameters:

	f <- k + p*k # weights -1 + means +1 => start of alpha
	f1 <- f # end of alpha if alpha uniform
	f1.1 <- f1 +1L #start of D. if alpha unif.
	f1.2 <- f1 + p # end of D. if D. uniform and alpha uniform
	f1.3 <- f1 + p*k # end D. if D. var and alpha unif.
	
	f2.1 <- f1 + k # start of D. if alpha varialbe
	f2.2 <- f2.1 + p # end of D. if D. uniform and alpha variable
	f2.3 <- f2.1 + p*k # end of D. if D.var and alpha var

	#only important ones are f1.2, f1.3, f2.2, f2.3


	model <- match.arg(model)

	par. <- switch(model,

	"EII" = {par.[f] <- log(par.[f])
		 par.},

	"VII" = {par.[f:f2.1] <- log(par.[f:f2.1])
		 par.},

	"EEI" = ,
	"EEE" = {par.[f:f1.1] <- log(par.[f:f1.1])
		 par.},

	"VEI" = ,
	"VEE" = {par.[f:f2.2] <- log(par.[f:f2.2])
		 par.},

	"EVI" = ,
	"EVV" = {par.[f:f1.3] <- log(par.[f:f1.3])
		 par.},

	"VVI" = ,
	"VVV" = {par.[f:f2.3] <- log(par.[f:f2.3])
		 par.}
	)

	par.
}
