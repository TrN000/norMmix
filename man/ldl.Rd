\name{ldl}
\alias{ldl}
\title{LDL' Cholesky Decomposition}
\description{
  Simple (but not too simple) \R implementation of the (square root free)
  \eqn{LDL'} Choleksy decomposition.
}
\usage{
ldl(m)
}
\arguments{
  \item{m}{positive semi-definite square matrix, say of dimension \eqn{n
      \times n}{n x n}.}
}
%% \details{
%% }
\value{
  a \code{\link{list}} with two components
 \item{L}{a lower triangular matrix with diagonal entries 1.}
 \item{D}{numeric vector, the \emph{diagonal}
   \eqn{d_{1,1},d_{2,2},\dots,d_{n,n}} of the diagonal matrix \eqn{D}.}
}
%% \references{
%% }
\seealso{
  \code{\link{chol}()} in base \R, or also a \dQuote{generalized LDL}
  decomposition, the Bunch-Kaufman, \code{\link[Matrix]{BunchKaufman}()}
  in (\sQuote{Recommended}) package \CRANpkg{Matrix}.
}
\examples{
(L <- rbind(c(1,0,0), c(3,1,0), c(-4,5,1)))
D <- c(4,1,9)
FF <- L \%*\% diag(D) \%*\% t(L)
FF
LL <- ldl(FF)
stopifnot(all.equal(L, LL$L),
          all.equal(D, LL$D))

## rank deficient :
FF0 <- L \%*\% diag(c(4,0,9)) \%*\% t(L)
((L0 <- ldl(FF0))) #  !! now fixed with the  if(Di == 0) test
## With the "trick", it works:
stopifnot(all.equal(FF0,
                    L0$L \%*\% diag(L0$D) \%*\% t(L0$L)))
## [hint: the LDL' is no longer unique when the matrix is singular]

system.time(for(i in 1:10000) ldl(FF) ) # ~ 0.2 sec

(L <- rbind(c( 1, 0, 0, 0),
            c( 3, 1, 0, 0),
            c(-4, 5, 1, 0),
            c(-2,20,-7, 1)))
D <- c(4,1, 9, 0.5)
F4 <- L \%*\% diag(D) \%*\% t(L)
F4
L4 <- ldl(F4)
stopifnot(all.equal(L, L4$L),
          all.equal(D, L4$D))

system.time(for(i in 1:10000) ldl(F4) )


## rank deficient :
F4.0 <- L \%*\% diag(c(4,1,9,0)) \%*\% t(L)
((L0 <- ldl(F4.0)))
stopifnot(all.equal(F4.0,
                    L0$L \%*\% diag(L0$D) \%*\% t(L0$L)))

F4_0 <- L \%*\% diag(c(4,1,0,9)) \%*\% t(L)
((L0 <- ldl(F4_0)))
stopifnot(all.equal(F4_0,
                    L0$L \%*\% diag(L0$D) \%*\% t(L0$L)))


## Large
mkLDL <- function(n, rF = function(n) sample.int(n), rFD = function(n) 1+ abs(rF(n))) {
    L <- diag(nrow=n)
    L[lower.tri(L)] <- rF(n*(n-1)/2)
    list(L = L, D = rFD(n))
}

(LD <- mkLDL(17))

chkLDL <- function(n, ..., verbose=FALSE, tol = 1e-14) {
    LD <- mkLDL(n, ...)
    if(verbose) cat(sprintf("n=\%3d ", n))
    n <- length(D <- LD$D)
    L <- LD$L
    M <- L \%*\% diag(D) \%*\% t(L)
    r <- ldl(M)
    stopifnot(exprs = {
        all.equal(M,
                  r$L \%*\% diag(r$D) \%*\% t(r$L), tol=tol)
        all.equal(L, r$L, tol=tol)
        all.equal(D, r$D, tol=tol)
    })
    if(verbose) cat("[ok]\n")
    invisible(list(LD = LD, M = M, ldl = r))
}

(chkLDL(7))

N <- 99 ## test  N  random cases
set.seed(101)
for(i in 1:N) {
    cat(sprintf("i=\%3d, ",i))
    chkLDL(rpois(1, lambda = 20), verbose=TRUE)
}

system.time(chkLDL( 500)) # 0.62

try( ## this almost never "works":
system.time(chkLDL( 500, rF = rnorm, rFD = function(n) 10 + runif(n))) # 0.64
)

if(interactive())
   system.time(chkLDL( 600)) # 1.09
## .. then it grows quickly for (on nb-mm4)
## for n = 1000  it typically *fails*: The matrix M  is typically very ill conditioned
## does not depend much on the RNG ?

"==> much better conditioned L and hence M : "
set.seed(120)
L <- as(Matrix::tril(toeplitz(exp(-(0:999)/50))), "matrix")
dimnames(L) <- NULL
D <- 10 + runif(nrow(L))
M <- L \%*\% diag(D) \%*\% t(L)
rcond(L) # 0.010006 !
rcond(M) # 9.4956e-5
if(FALSE) # ~ 4-5 sec
   system.time(r <- ldl(M))
}
