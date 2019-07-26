#### functions handling parameter manipulation. par2nM nM2par etc

### for info on models see Celeux and Govaert 1995

#' @include Cholesky.R
#' @include weight.R


## map lower.tri to vec
ld. <- function(mat){
	x <- mat[lower.tri(mat, diag=FALSE)]
}

## map vec to lower.tri
dl. <- function(d,x,p){
	mat <- diag(1,p)
	mat[lower.tri(mat,diag=FALSE)] <- x
	mat %*% diag(d) %*% t(mat)
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
#' @seealso n2p
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

	stopifnot( isTRUE(all.equal(sum(w),1)), (length(w)==k), 
		  is.numeric(w), is.finite(w) )


	# mu

	stopifnot( (dim(mu)==c(p,k)), (is.numeric(mu)),
		  is.matrix(mu), is.finite(mu) )


	# Sigma

	stopifnot( (dim(sig)==c(p,p,k)), (is.numeric(sig)),
		  length(dim(sig))==3, is.array(sig) )
	stopifnot( isTRUE( all( apply(sig,3, function(j) (ldl(j)$Diag >= 0 )))))


	#output vector of parameter values

	c(
	  w <- switch(trafo, #weights either logit or centered log ratio
		"logit" = logit(w),

    		"clr1" = clr1(w),

		stop("Error in weight trafo, ",trafo)
	    	),
	  mu, #means
	  Sigma <- switch(model, #model dependent covariance values
		"EII" = log(sig[1,1,1]),

		"VII" = log(sig[1,1,]),

		"EEI" = {D.temp <- diag(sig[,,1])
			alpha <- log(prod(D.temp)^(1/p))
			D. <- log(D.temp) - alpha
			c(alpha, D.)},

		"VEI" = {alpha <- (apply(sig,3,function(j) prod(diag(j))^(1/p)))
			D. <- (diag(sig[,,1]))/alpha[1]
			c(log(alpha), log(D.))},
			  
		"EVI" = {alpha <- log(prod(diag(sig[,,1]))^(1/p))
			D. <- log(apply(sig,3,diag)) - alpha
			c(alpha,D.)},

		"VVI" = {alpha <- apply(sig,3, function(j) det(j)^(1/p)) 
			D.temp <- apply(sig,3,diag)
			D. <- D.temp %*% diag(1/alpha) # this is fastest?? https://stackoverflow.com/questions/20596433/how-to-divide-each-row-of-a-matrix-by-elements-of-a-vector-in-r
			c(log(alpha),log(D.))},

		"EEE" = {alpha <- prod( ldl(sig[,,1])$Diag )^(1/p)
			A. <- ldl(sig[,,1])
			c(log(alpha), log(A.$Diag/alpha), ld.(A.$L))},

		"VEE" = {alpha <- apply(sig,3, function(j) prod(ldl(j)$Diag )^(1/p))
			A. <- ldl(sig[,,1])
			c(log(alpha), log(A.$Diag/alpha[1]), ld.(A.$L))},

		"EVV" = {alpha <- prod( ldl(sig[,,1])$Diag )^(1/p)
			D. <- apply(sig,3, function(j) ldl(j)$Diag/alpha)
			L. <- apply(sig,3, function(j) ld.( ldl(j)$L ))
			c(log(alpha), log(D.), L.)},

		"VVV" = {alpha <- apply(sig,3, function(j) prod(ldl(j)$Diag )^(1/p))
			D.temp <- apply(sig,3, function(j) ldl(j)$Diag)
			D. <- D.temp %*% diag(1/alpha)
			L. <- apply(sig,3, function(j) ld.( ldl(j)$L ))
			c(log(alpha), log(D.), L.)},

		stop("invalid argument in 'model'")
		)
	)
}


#' wrapper function for nMm objs in zmarrwandMm
#'
#' \code{n2p} returns same as nMm2par with clr1
#'
#' @export

n2p <- function(obj) nMm2par(obj, trafo="clr1", model=obj$model)




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

	p <- as.integer(p)
	k <- as.integer(k)


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

	#only important ones are f1.2, f1.3, f2.2, f2.3


	par. <- switch(model,

	"EII" = {par.[f] <- exp(par.[f])
		 par.},

	"VII" = {par.[f:f2]<- exp(par.[f:f2])
		 par.},

	"EEI" = ,
	"EEE" = {par.[f:f11] <- exp(par.[f:f11])
		 par.},

	"VEI" = ,
	"VEE" = {par.[f:f21] <- exp(par.[f:f21])
		 par.},

	"EVI" = ,
	"EVV" = {par.[f:f12] <- exp(par.[f:f12])
		 par.},

	"VVI" = ,
	"VVV" = {par.[f:f22] <- exp(par.[f:f22])
		 par.},

	stop("Error in exp switch statement in par2nMm")
	)

	if (k ==1) {w.temp <- vector()}
	else {w.temp <- par.[1:(k-1)]}
	w <- switch(trafo,
	"logit" = logitinv(w.temp),

	"clr1" = clr1inv(w.temp), 
	)

	mu <- matrix(par.[k:(k+p*k-1)], p, k)



	Sigma <- switch(model,

	# diagonal cases
	"EII" = {lambda <- par.[f]

		 sig <- array( rep(diag(lambda, p),k), c(p,p,k) )},

	"VII" = {lambda <- par.[f:f2]

		 sig <- array(unlist(lapply( lambda, function(j) diag(j,p) )), c(p,p,k)) },

	"EEI" = {lambda <- par.[f]

		 D.temp <- par.[f1.1:f11]

		 sig <- array( rep(diag(lambda*D.temp),k), c(p,p,k) )},

	"VEI" = {lambda <- par.[f:f2]

		 D.temp <- par.[f2.1:f21]

		 D. <- tcrossprod(D.temp,lambda)

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	"EVI" = {lambda <- par.[f]

		 D. <- matrix(par.[f1.1:f12],p,k)*lambda

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	"VVI" = {lambda <- par.[f:f2]

		 D.temp <- matrix(par.[f2.1:f22],p,k)

		 D. <- D.temp %*% diag(lambda)

		 sig <- array( apply(D.,2, diag), c(p,p,k)) },

	# variable cases

	"EEE" = {lambda <- par.[f]
		 D.temp <- par.[f1.1:f11]
		 D. <- D.temp*lambda

		 L. <- par.[f11.1:f111]
		 A. <- dl.(D.,L.,p)

		 sig <- array(0, c(p,p,k))
		 for (i in 1:k){
			 sig[,,i] <- A.
		 }
		 sig},

	"VEE" = {lambda <- par.[f:f2]
		 D.temp <- par.[f2.1:f21]

		 D. <- tcrossprod(D.temp,lambda)

		 f3 <- (p*(p-1)/2)

		 L. <- par.[f21.1:f211]

		 sig <- array(0, c(p,p,k))
		 for (i in 1:k){
			 sig[,,i] <- dl.(D.[,i],L.,p)
		 }
		 sig},

	"EVV" = {lambda <- par.[f]

		 D. <- matrix(par.[f1.1:f12],p,k) * lambda

		 f3 <- (p*(p-1)/2)

		 L.temp <- matrix(par.[f12.1:f121],f3,k)

		 sig <- array(0, c(p,p,k))
		 for (i in 1:k) {
			 sig[,,i] <- dl.(D.[,i],L.temp[,i],p)
		 }
		 sig},

	"VVV" = {lambda <- par.[f:f2]

		 D.temp <- matrix(par.[f2.1:f22],p,k)

		 D. <- D.temp %*% diag(lambda)

		 f3 <- (p*(p-1)/2)

		 L.temp <- matrix(par.[f22.1:f221],f3,k)

		 sig <- array(0, c(p,p,k))
		 for (i in 1:k) {
			 sig[,,i] <- dl.(D.[,i],L.temp[,i],p)
		 }
		 sig},
	stop("error in Sigma switch statement")
	)


	name <- sprintf("model = %s , clusters = %s", model, k)


	structure(
		  name = name,
		  class = "norMmix",
		  list( mu=mu, Sigma=Sigma, weight=w, k=k, dim=p , model=model)
		  )

}



