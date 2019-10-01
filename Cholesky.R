### LDL Cholesky decomposition, not particularly efficient, but works for now


##' MM:  Now also "works" in the rank deficient case!
ldl <- function(m) {
    stopifnot(is.matrix(m), is.numeric(m),
              (n <- nrow(m)) == ncol(m), n >= 1)
    D <- Zn <- numeric(n)# = rep(0,n)
    L <- m
    for (i in 1:n) {
        D[i] <- Di <- (mi <- m[,i])[i]
        if(Di != 0) {
            L[,i] <- Li <- mi/Di
            m <- m - Di * tcrossprod(Li)
        } else  ## Di == 0:
            L[,i] <- 0

        m[i,] <- m[,i] <- Zn
    }
    list(L=L, D=D)
}


# inverse to ldl

ldlinv <- function(D,L) {
    stopifnot(is.matrix(L), is.vector(D), identical(dim(L), c(p, p <- length(D))))
    L %*% diag(D, nrow=p) %*% t(L)
}

## ==> Examples at end of  ../man/ldl.Rd
##                         -------------
