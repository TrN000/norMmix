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


#(L <- rbind(c(1,0,0), c(3,1,0), c(-4,5,1)))
#D <- c(4,1,9)
#FF <- L %*% diag(D) %*% t(L)
#FF
#LL <- ldl(FF)
#stopifnot(all.equal(L, LL$L),
#          all.equal(D, LL$D))
#
### rank deficient :
#FF0 <- L %*% diag(c(4,0,9)) %*% t(L)
#((L0 <- ldl(FF0))) #  !! now fixed with the  if(Di == 0) test
### With the "trick", it works:
#stopifnot(all.equal(FF0,
#                    L0$L %*% diag(L0$D) %*% t(L0$L)))
### [hint: the LDL' is no longer unique when the matrix is singular]
#
#system.time(for(i in 1:10000) ldl(FF) )
###  user  system elapsed
### 0.174   0.031   0.205 < LT
### 0.161   0.000   0.162 < MM
#
#(L <- rbind(c( 1, 0, 0, 0),
#            c( 3, 1, 0, 0),
#            c(-4, 5, 1, 0),
#            c(-2,20,-7, 1)))
#D <- c(4,1, 9, 0.5)
#F4 <- L %*% diag(D) %*% t(L)
#F4
#L4 <- ldl(F4)
#stopifnot(all.equal(L, L4$L),
#          all.equal(D, L4$D))
#
#system.time(for(i in 1:10000) ldl(F4) )
#
#
### rank deficient :
#F4.0 <- L %*% diag(c(4,1,9,0)) %*% t(L)
#((L0 <- ldl(F4.0)))
#stopifnot(all.equal(F4.0,
#                    L0$L %*% diag(L0$D) %*% t(L0$L)))
#
#F4_0 <- L %*% diag(c(4,1,0,9)) %*% t(L)
#((L0 <- ldl(F4_0)))
#stopifnot(all.equal(F4_0,
#                    L0$L %*% diag(L0$D) %*% t(L0$L)))
#
#
### Large
#mkLDL <- function(n, rF = function(n) sample.int(n), rFD = function(n) 1+ abs(rF(n))) {
#    L <- diag(nrow=n)
#    L[lower.tri(L)] <- rF(n*(n-1)/2)
#    list(L = L, D = rFD(n))
#}
#
#(LD <- mkLDL(17))
#
#chkLDL <- function(n, ..., verbose=FALSE, tol = 1e-14) {
#    LD <- mkLDL(n, ...)
#    if(verbose) cat(sprintf("n=%3d ", n))
#    n <- length(D <- LD$D)
#    L <- LD$L
#    M <- L %*% diag(D) %*% t(L)
#    r <- ldl(M)
#    stopifnot(exprs = {
#        all.equal(M,
#                  r$L %*% diag(r$D) %*% t(r$L), tol=tol)
#        all.equal(L, r$L, tol=tol)
#        all.equal(D, r$D, tol=tol)
#    })
#    if(verbose) cat("[ok]\n")
#    invisible(list(LD = LD, M = M, ldl = r))
#}
#
#(chkLDL(7))
#
### test a 1000 random cases
#set.seed(101)
#for(i in 1:1000) {
#    cat(sprintf("i=%3d, ",i))
#    chkLDL(rpois(1, lambda = 20), verbose=TRUE)
#}
#
#system.time(chkLDL( 500)) # 0.62
#
#try( ## this almost never "works":
#system.time(chkLDL( 500, rF = rnorm, rFD = function(n) 10 + runif(n))) # 0.64
#)
#
#system.time(chkLDL( 600)) # 1.09
### .. then it grows quickly for (on nb-mm4)
### for n = 1000  it typically *fails*: The matrix M  is typically very ill conditioned
### does not depend much on the RNG ?
#
#"==> much better conditioned L and hence M : "
#L <- as(tril(toeplitz(exp(-(0:999)/50))), "matrix")
#dimnames(L) <- NULL
#D <- 10 + runif(nrow(L))
#M <- L %*% diag(D) %*% t(L)
#rcond(L) # 0.01000 !!
#rcond(M) # 9.478e-5 .
#system.time(r <- ldl(M))
###  user  system elapsed
### 5.207   0.153   5.412
#
#
#
## with isSymmetric
## system.time(  {for (i in 1:10000){ldl(FF)}}  )
##    user  system elapsed
##   2.236   0.000   2.235
#
#
## system.time(  {for (i in 1:1000000){ldl(FF)}}  )
##    user  system elapsed
##  24.404   0.008  24.434
#
#
