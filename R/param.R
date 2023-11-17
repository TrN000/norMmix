#### functions handling parameter manipulation. par2nM nM2par etc

### for information on models, see Celeux and Govaert 1995


## map lower.tri to vec
ld. <- function(mat) mat[lower.tri(mat, diag=FALSE)]

## map vec to lower.tri
dl. <- function(d,x,p) {
    mat <- diag(1,p)
    mat[lower.tri(mat,diag=FALSE)] <- x
    mat %*% diag(d) %*% t(mat) ## << FIXME: make *fast* for larger 'mat'
}


## centered log ratio
clr1 <- function(w) {
    stopifnot(length(w) >= 1L, w >= 0, all.equal(sum(w), 1))
    # calculate clr1
    ln <- log(w)
    ## log(0) == -Inf "kills"  mean(ln) etc:
    if(any(w0 <- !w)) ln[w0] <- -709 # < log(.Machine$double.xmin) = -708.3964
    ## return:
    ln[-1L] - mean(ln)
}

.logMax <- log(2) * .Machine$double.max.exp # = 709.7827 = log(.Machine$double.xmax)

clr1inv <- function(lp) {
    ## stopifnot(is.numeric(lp))
    if(!length(lp)) return(1)
    ## calc weights
    lp <- c(-sum(lp), lp) # = (lp_1,..., lp_m)
    if((mlp <- max(lp))*length(lp) < .logMax) { ## normal case
        p1 <- exp(lp)
    } else { ## mlp = max(lp) >= .logMax, where exp(lp) would overflow
        p1 <- exp(lp - mlp)
    }
    p1/sum(p1)
}

norModels <- eval(formals(norMmix)$model)
## ==  c("EII","VII","EEI","VEI","EVI",
##       "VVI","EEE","VEE","EVV","VVV")

# nMm2par returns vector of parameters of norMmix objects
#
# This transformation forms a vector from the parameters of a normal
# mixture. These consist of weights, means and covariance matrices.
# Cov mats are given as D and L from the LDLt decomposition
#
# see also n2p
#
# obj:   list containing sig= covariance matrix array, mu= mean vector matrix,
#        w= weights, k= number of components, p= dimension
# model: one of "EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVV" or "VVV"
nMm2par <- function(obj
                  , model = c("EII","VII","EEI","VEI","EVI",
                              "VVI","EEE","VEE","EVV","VVV")
                  , meanFUN = mean.default
                  , checkX = FALSE
                    ) {
    w   <- obj$weight
    mu  <- obj$mu
    sig <- obj$Sigma
    p <- obj$dim
    k <- obj$k

    model <- match.arg(model)

    av <- match.fun(meanFUN)

    ## Checks -----------------
    # weights:
    stopifnot(is.numeric(w), all.equal(sum(w),1), length(w) == k)
    # mu:
    stopifnot(is.matrix(mu), dim(mu) == c(p,k), is.numeric(mu), is.finite(mu))
    # Sigma: -- FIXME: once we allow *compact* sig ==> the checks get partly get much cheaper
    stopifnot(is.array(sig), is.numeric(sig), length(dim(sig)) == 3L
            , dim(sig) == c(p,p,k)
            , apply(sig, 3L,
                    if(checkX) ## Expensive(!)
                         function(j) ldl(j)$D >= 0
                    else function(j) diag(j)  >= 0)
              )

    ## output vector of parameter values

    c(# weights 'w' (= \pi_j ):
      clr1(w),
      mu, # means
      ## Sigma :
      switch(model, # model dependent covariance values
        "EII" = {
            lD. <- log(apply(sig,3, function(j) ldl(j)$D))
            av(lD.)
            },

        "VII" = {
            lD. <- log(apply(sig,3, function(j) ldl(j)$D))
            apply(lD., 2, av) # faster for mean: colMeans(lD.)
            },

        "EEI" = {
            D. <- apply(sig,3, function(j) ldl(j)$D)
            D. <- apply(D.,1, function(j) av(log(j)))
            alpha <- mean(D.)
            D. <- D.-alpha
            c(alpha, D.[-1])
            },

        "VEI" = {
            D. <- apply(sig,3, function(j) ldl(j)$D)
            alpha <- apply(D.,2, function(j) av(log(j)))
            D.    <- apply(D.,1, function(j) av(log(j)))
            D. <- D.-mean(D.)
            c(alpha, D.[-1])
            },

        "EVI" = {
            D. <- apply(sig,3, function(j) ldl(j)$D)
            alpha <- av(log(D.))
            D. <- log(D.)
            D. <- apply(D.,2, function(j) j-av(j))
            c(alpha,D.[-1,])
            },

        "VVI" = {
            D. <- apply(sig,3, function(j) log(ldl(j)$D))
            alpha <- apply(D.,2, av) # faster for mean: colMeans(lD.)
            D. <- apply(D.,2, function(j) j-av(j))
            c(alpha, D.[-1,])
            },

        "EEE" = {
            D. <- apply(sig,3, function(j) log(ldl(j)$D))
            D. <- apply(D.,1, av)  # faster for mean: rowMeans(lD.) (?)
            alpha <- av(D.)
            D. <- D.-av(D.)
            L. <- apply(sig,3, function(j) ld.(ldl(j)$L))
            L. <- matrix(L., p*(p-1)/2, k)
            L. <- apply(L.,1, av)
            c(alpha, D.[-1], L.)
            },

        "VEE" = {
            D. <- apply(sig,3, function(j) log(ldl(j)$D) )
            alpha <- apply(D.,2, av)
            D. <- apply(D.,1, av)
            D. <- D.-av(D.)
            L. <- apply(sig,3, function(j) ld.(ldl(j)$L))
            L. <- matrix(L., p*(p-1)/2, k)
            L. <- apply(L.,1, av)
            c(alpha, D.[-1], L.)
            },

        "EVV" = {
            D. <- apply(sig,3, function(j) log(ldl(j)$D))
            alpha <- av(D.)
            D. <- apply(D.,2, function(j) j-av(j))
            L. <- apply(sig,3, function(j) ld.(ldl(j)$L))
            c(alpha, D.[-1,], L.)
            },

        "VVV" = {
            D. <- apply(sig,3, function(j) log(ldl(j)$D))
            alpha <- apply(D.,2, av)
            D. <- apply(D.,2, function(j) j-av(j))
            L. <- apply(sig,3, function(j) ld.(ldl(j)$L))
            c(alpha, D.[-1,], L.)
            },

        stop("invalid 'model': ", model)
        )
    )
}


# transform of parameter vector to normal mixture
#
# par2nMm returns list containing weight, mu, Sigma, k, dim
#
# this is the inverse function to nMm2par. Given a numeric vector
# dimension and cluster number this function reconstructs a normal mixture
# object.
#
# par:   numeric vector of parameters (alpha, D., L.)
# p:     dimension of space
# model: See description
#
# returns list: list(weight=w, mu=mu, Sigma=Sigma, k=k, dim=p)
par2nMm <- function(par, p, k
                  , model = c("EII","VII","EEI","VEI","EVI",
                              "VVI","EEE","VEE","EVV","VVV")
                  , name = sprintf("model = %s , components = %s", model, k)
                    ) {
    model <- match.arg(model)

    p <- as.integer(p)
    k <- as.integer(k)

    # start of relevant parameters:

    f <- k + p*k # weights -1 + means +1 => start of alpha

    f1 <- f      # end of alpha if uniform
    f2 <- f+k-1L # end of alpha if var

    f1.1 <- f1 +1L # start of D. if alpha unif.
    f2.1 <- f1 + k # start of D. if alpha variable

    p1 <- p - 1L
    f11 <- f1 + p1   # end of D. if D. uniform and alpha uniform
    f12 <- f1 + p1*k # end of D. if D. var     and alpha unif.
    f21 <- f2 + p1   # end of D. if D. uniform and alpha variable
    f22 <- f2 + p1*k # end of D. if D. var     and alpha var

    f11.1 <- f11 +1L # start of L if alpha unif D unif
    f21.1 <- f21 +1L # start of L if alpha var  D unif
    f12.1 <- f12 +1L # start of L if alpha unif D var
    f22.1 <- f22 +1L # start of L if alpha var  D var

    pC2 <- (p*p1) %/% 2L # == choose(p, 2) == p(p-1)/2
    f111 <- f11 +   pC2 # end of L if alpha unif D unif
    f211 <- f21 +   pC2 # end of L if alpha var  D unif
    f121 <- f12 + k*pC2 # end of L if alpha unif D var
    f221 <- f22 + k*pC2 # end of L if alpha var  D var

    # only important ones are f1.2, f1.3, f2.2, f2.3

    w <- clr1inv(par[seq_len(k-1L)])

    mu <- matrix(par[k:(f-1L)], p, k)

### FIXME: Alternatively, instead of Sigma, compute  chol(Sigma) = D^{1/2} L'  as Sigma = LDL'

    Sigma <- switch(model,
    # diagonal cases
    "EII" = {
        alpha <- exp(par[f])
        array( rep(diag(alpha, p),k), c(p,p,k) )
        },

    "VII" = {
        alpha <- exp(par[f:f2])
        array(unlist(lapply( alpha, function(j) diag(j,p) )), c(p,p,k))
        },

    "EEI" = {
        alpha <- par[f]
        D. <- par[f1.1:f11]
        D. <- c(-sum(D.), D.)
        D. <- D.-mean(D.)
        array( rep(diag(exp(alpha+D.)),k), c(p,p,k) )
        },

    "VEI" = {
        alpha <- par[f:f2]
        D. <- par[f2.1:f21]
        D. <- c(-sum(D.), D.)
        D. <- matrix(D.+rep(alpha, each=p), p, k)
        D. <- exp(D.)
        array( apply(D.,2, diag), c(p,p,k))
        },

    "EVI" = {
        alpha <- par[f]
        D. <- matrix(par[f1.1:f12],p-1,k)
        D. <- apply(D., 2, function(j) c(-sum(j), j))
        D. <- exp(D.+alpha)
        array( apply(D.,2, diag), c(p,p,k))
        },

    "VVI" = {
        alpha <- par[f:f2]
        D. <- matrix(par[f2.1:f22],p-1,k)
        D. <- apply(D., 2, function(j) c(-sum(j), j))
        D. <- exp(D.+rep(alpha, each=p))
        array(apply(D.,2, diag), c(p,p,k))
    },

    # variable cases

    "EEE" = {
        alpha <- par[f]
        D. <- par[f1.1:f11]
        D. <- c(-sum(D.), D.)
        D. <- exp(D.+alpha)
        L. <- par[f11.1:f111]
        A. <- dl.(D.,L.,p)
        ## sig :=
        array(rep(A., times=k), c(p,p,k))
    },

    "VEE" = {
        alpha <- par[f:f2]
        D. <- par[f2.1:f21]
        D. <- c(-sum(D.), D.)
        D. <- exp(matrix(D.+rep(alpha, each=p), p, k))
        L. <- par[f21.1:f211]
        sig <- array(0, c(p,p,k))
        for (i in 1:k){
            sig[,,i] <- dl.(D.[,i],L.,p)
        }
        sig
    },

    "EVV" = {
        #par[f:f12] <- exp(par[f:f12])
        alpha <- par[f]
        D. <- matrix(par[f1.1:f12],p-1,k)
        D. <- apply(D., 2, function(j) c(-sum(j), j))
        D. <- exp(D.+alpha)
        L.temp <- matrix(par[f12.1:f121], pC2,k)
        sig <- array(0, c(p,p,k))
        for (i in 1:k) {
            sig[,,i] <- dl.(D.[,i],L.temp[,i],p)
        }
        sig
    },

    "VVV" = {
        #par[f:f22] <- exp(par[f:f22])
        alpha <- par[f:f2]
        D. <- matrix(par[f2.1:f22],p-1,k)
        D. <- apply(D., 2, function(j) c(-sum(j), j))
        D. <- exp(D.+rep(alpha,each=p))
        L.temp <- matrix(par[f22.1:f221], pC2,k)
        sig <- array(0, c(p,p,k))
        for (i in 1:k) {
            sig[,,i] <- dl.(D.[,i],L.temp[,i],p)
        }
        sig
    },
    stop("error in Sigma switch statement")
    )

    structure(
        name = name,
        class = "norMmix",
        list( mu=mu, Sigma=Sigma, weight=w, k=k, dim=p, model=model)
        )
}


nparSigma <- function(k, p,
                      model=c("EII","VII","EEI","VEI","EVI",
                              "VVI","EEE","VEE","EVV","VVV"))
    switch(match.arg(model),
        "EII" = 1L,
        "VII" = k,
        "EEI" = 1L+   (p-1L),
        "VEI" = k +   (p-1L),
        "EVI" = 1L+ k*(p-1L),
        "VVI" = k + k*(p-1L),
        "EEE" = 1L+   (p-1L) +   (p*(p-1L)) %/% 2L,
        "VEE" = k +   (p-1L) +   (p*(p-1L)) %/% 2L,
        "EVV" = 1L+ k*(p-1L) + k*(p*(p-1L)) %/% 2L,
        "VVV" = k + k*(p-1L) + k*(p*(p-1L)) %/% 2L # == k * p*(p+1)/2
        )


dfnMm <- function(k, p,
                   model=c("EII","VII","EEI","VEI","EVI",
                           "VVI","EEE","VEE","EVV","VVV")
                   ) {
    stopifnot(is.finite(k <- as.integer(k)), is.finite(p <- as.integer(p)))
    model <- match.arg(model)
    w <- k-1L
    mu <- p*k
    sig <- nparSigma(k, p, model)
    ## return
    w + mu + sig
}


npar <- function(object, ...) {UseMethod("npar")}


parcond <- function(x,
                    k,
                    model=c("EII","VII","EEI","VEI","EVI",
                            "VVI","EEE","VEE","EVV","VVV")
                    ) {
    stopifnot(is.matrix(x), length(k) == 1L, k == as.integer(k), k >= 1)
    n <- nrow(x)
    p <- ncol(x)
    model <- match.arg(model)
    n / dfnMm(k, p, model=model)
}

