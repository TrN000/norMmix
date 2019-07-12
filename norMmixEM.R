### EM-algorithm for norMmix datasets

## following the description in McLachlan & Peel (2000) p.82

## nomenclature used throughout:
## n sample size
## k number of clusters
## p dimension of sample

library(mixtools)


initPar <- function(tau, mu, Sigma, weight){

	#checking inputs here

	par <- structure(class="par", .Data=list(tau=tau, mu=mu, Sigma=Sigma, weight=weightr))

}


# M step

mstep.eMm <- function( y, tau, mu, Sigma, weight ){
	
	if (!is.matrix(y)) stop("in mstep y is not matrix")

	#  calculating T_i

	T1 <- rowSums(tau) # vec of length k

	T2 <- tau %*% y # matrix of size k x p

	T3 <- array(0, dim=c(p,p,k))

	for (i in 1:k){
		for (j in 1:n){
		T3[,,i] <- T3[,,i] + (t(tau[i,]* y[j,])%*%y[j,])
		}
	}


	#mu

	mu[i,] <- T2/T1

	#sigma

	Sigma <- array(0, c(p,p,k))

	for (i in 1:k){
		Sigma[,,i] <- ( T3[,,i]-T1[i]^-1 *(t(T2[,i])%*%T2[,i]) ) / T1[i]
	}

	# weight

	weight <- T1/ncol(y)

	##return params

	list(tau, weight=weight, mu=mu, Sigma=Sigma)

}


estep.nMm <- function(y, par){

	mu <- par$mu
	Sigma <- par$Sigma
	weight <- par$weight

	phi <- matrix(0,k,n)
	
	for (i in 1:k){
	phi[i,] <- dmvnorm(y, mu=mu[i,], sigma=Sigma[,,i]) # a k x n vec
	}

	T5 <- rep(0,n) 
	T5 <- colSums( phi %*% diag(weight))^-1

	tau <- (weight %*% phi) %*% T5
	

	## return params
	par$tau <- tau
}

	
