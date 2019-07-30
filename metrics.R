### implements several measures of distribution differences to be used by 
### metric.norMmix(), which in turn is used to test norMmixMLE


#' @include norMmix.R


hell <- function(nMm1, nMm2) {

	is.norMmix(nMm1)
	is.norMmix(nMm2)

	k1 <- nMm1$k
	k2 <- nMm2$k

	stopifnot(k1==k2)

	p1 <- nMm1$dim
	p2 <- nMm2$dim

	mu1 <- nMm1$mu
	mu2 <- nMm2$mu

	sig1 <- nMm1$Sigma
	sig2 <- nMm2$Sigma

	order. <- match.norMmix(nMm1, nMm2)

	res <- vector()

	for (i in 1:k) {
	res[i] <- 1- det(sig1[,,i])^(1/4)*det(sig2[order.[i]])^(1/4)/sqrt(det((sig1[,,i]*sig2[,,order.[i]])/2))*
		exp(t((mu1[,i]-mu[,order.[i]]))%*%((sig1[,,i]+sig2[,,order.[i]]))%*%(mu1[,i]-mu[,order.[i]])/8)
	}

	return(prod(res))
}
