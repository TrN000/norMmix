#### EM-algorithm for norMmix datasets

## following the description in McLachlan & Peel (2000) p.82
## M step
mstep.nMm <- function(x, tau) {
    stopifnot(is.matrix(x), is.matrix(tau), nrow(tau) == (n <- nrow(x)))
    # x   is  n x p matrix
    p <- ncol(x)
    k <- ncol(tau) # tau is  n x k

    ## calculating T_i
    T1 <- colSums(tau) # vec of length k, integer sum, stable
    T2 <- t(tau) %*% x # matrix of size k x p
    T3 <- array(0, dim=c(p,p,k))
    for (i in 1:k){
        for (j in 1:n){
        T3[,,i] <- T3[,,i] + tcrossprod((tau[j,i]* x[j,]),x[j,]) # could improve, since we know tau={0,1}^n*k
        }
    }

    mu <- t(T2)/rep(T1,each=p)
    Sigma <- array(0, c(p,p,k))
    for (i in 1:k){
        Sigma[,,i] <- ( T3[,,i] - tcrossprod(T2[i,],T2[i,])/T1[i]) / T1[i]
    }
    weight <- T1/nrow(x)

    ##return params
    norMmix(mu=mu, Sigma=Sigma, weight=weight, model="VVV")
}
