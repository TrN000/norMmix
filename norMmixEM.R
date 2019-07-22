### EM-algorithm for norMmix datasets

## following the description in McLachlan & Peel (2000) p.82

## nomenclature used throughout:
## n sample size
## k number of clusters
## p dimension of sample

#library(mixtools)


#initPar <- function(tau, mu, Sigma, weight){
#
#	#checking inputs here
#
#	par <- structure(class="par", .Data=list(tau=tau, mu=mu, Sigma=Sigma, weight=weightr))
#
#}
#

# M step

mstep.nMm <- function( x, tau, mu, Sigma, weighti, k ,p){
	
	if (!is.matrix(x)) stop("in mstep x is not matrix")

	n <- nrow(x)

	# x is n x p matrix
	# tau is n x k

	#  calculating T_i

	T1 <- colSums(tau) # vec of length k

	T2 <- t(tau) %*% x # matrix of size k x p

	T3 <- array(0, dim=c(p,p,k))

	for (i in 1:k){
		for (j in 1:n){
		T3[,,i] <- T3[,,i] + tcrossprod((tau[j,i]* x[j,]),x[j,])
		}
	}


	#mu

	mu <- t(T2)/T1

	#sigma

	Sigma <- array(0, c(p,p,k))

	for (i in 1:k){
		Sigma[,,i] <- ( T3[,,i]-T1[i]^(-1)*tcrossprod(T2[i,],T2[i,]) ) / T1[i]
	}

	# weight

	weight <- T1/nrow(x)

	##return params

	list(tau=tau, w=weight, mu=mu, Sigma=Sigma, k=k, dim=p)

}

# estep doesn't work yet, also not necessary.

#estep.nMm <- function(x, par){
#
#	mu <- par$mu
#	Sigma <- par$Sigma
#	weight <- par$weight
#	k <- par$k
#
#	phi <- matrix(0,k,n)
#	
#	for (i in 1:k){
#	phi[i,] <- dmvnorm(x, mu=mu[i,], sigma=Sigma[,,i]) # a k x n vec
#	}
#
#	T5 <- rep(0,n) 
#	T5 <- colSums( phi %*% diag(weight))^-1
#
#	tau <- (weight %*% phi) %*% T5
#	
#
#	## return params
#	par$tau <- tau
#	par
#}
#
	
