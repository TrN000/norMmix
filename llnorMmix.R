#### the llnorMmix function, calculating log likelihood for a given
#### parameter vector

## Author: Nicolas Trutmann 2019-07-06

#' @include param.R




{}



#' Log-likelihood of parameter vector given data
#'
#' \code{llnorMmix} returns scalar log-likelihood
#'
#'
#' description
#'
#' @param par. parameter vector
#' @param x sample matrix
#' @param k number of clusters
#' @param trafo either centered log ratio or logit
#' @param model assumed distribution model of normal mixture
#'
#' @export

llnorMmix <- function(par., x, k,
              trafo=c("clr1", "logit"),
              model=c("EII","VII","EEI","VEI","EVI",
                      "VVI","EEE","VEE","EVV","VVV")
              )
{
    stopifnot(is.matrix(x),
              length(k <- as.integer(k)) == 1, k >= 1)
    p <- ncol(x)
    x <- t(x) ## then only needed in   (x-mu[,i])^2  i=1..k

    # 2. transform

    trafo <- match.arg(trafo)
    model <- match.arg(model)

    l2pi <- log(2*pi)

    # 3. calc log-lik

    # get w

    w <- switch(trafo,
            "clr1" = {
                 if (k==1) 1
                     else clr1inv(par.[1:(k-1)])
                     },

            "logit" = {
                  if (k==1) 1
                      else logitinv(par.[1:(k-1)])
                      },

            stop("error in w switch in llnorMmix")
            )

#    if (!(sum(w)==1)) return(-Inf) ## FIXME allow tolerance

    # start of relevant parameters:

    f <- k + p*k # weights -1 + means +1 => start of alpha

    # get mu
    mu <- matrix(par.[k:(f-1L)], p,k)

    f1 <- f # end of alpha if uniform
    f2 <- f+k-1L # end of alpha if var

    f1.1 <- f1 +1L #start of D. if alpha unif.
    f2.1 <- f1 + k # start of D. if alpha varialbe

    f11 <- f1 + p -1# end of D. if D. uniform and alpha uniform
    f12 <- f1 + p*k -k# end D. if D. var and alpha unif.
    f21 <- f2 + p -1# end of D. if D. uniform and alpha variable
    f22 <- f2 + p*k -k# end of D. if D.var and alpha var

    f11.1 <- f11 +1L # start of L if alpha unif D unif
    f21.1 <- f21 +1L # start of L if alpha var D unif
    f12.1 <- f12 +1L # start of L if alpha unif D var
    f22.1 <- f22 +1L # start of L if alpha var D var

    f111 <- f11 + p*(p-1)/2 # end of L if alpha unif D unif
    f211 <- f21 + p*(p-1)/2 # end of L if alpha var D unif
    f121 <- f12 + k*p*(p-1)/2 # end of L if alpha unif D var
    f221 <- f22 + k*p*(p-1)/2 # end of L if alpha var D var


    # initialize f(x_i) i=1..n  vector of density values
    invl <- 0

    # calculate log-lik, see first case for explanation
    switch(model,
    "EII" = {
        alpha <- par.[f]
        invalpha <- exp(-alpha)# = 1/exp(alpha)
        for (i in 1:k) {
            rss <- colSums(invalpha*(x-mu[,i])^2)
            # this is vector of length n=sample size
            # calculates (x-mu)t * Sigma^-1 * (x-mu) for diagonal
            # cases.
            invl <- invl+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
            # adds likelihood of one component to invl
            # the formula in exp() is the log of likelihood
            # still of length n
        }
    },
    # hereafter differences are difference in dimension in alpha and D.
    # alpha / alpha[i] and D. / D.[,i]

    "VII" = {
        alpha <- par.[f:f2]
        for (i in 1:k) {
            rss <- colSums((1/exp(alpha[i]))*(x-mu[,i])^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
        }
    },

    "EEI" = {
        alpha <- par.[f]
        D. <- par.[f1.1:f11]
        D. <- c(-sum(D.),D.)
        D. <- D.-sum(D.)/p
        invD <- (1/exp(alpha+D.))
        for (i in 1:k) {
            rss <- colSums(invD*(x-mu[,i])^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
        }
    },

    "VEI" = {
        alpha <- par.[f:f2]
        D. <- par.[f2.1:f21]
        D. <- c(-sum(D.), D.)
        D. <- D.-sum(D.)/p
        for (i in 1:k) {
            rss <- colSums((1/exp(alpha[i]+D.))*(x-mu[,i])^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
        }
    },

    "EVI" = {
        alpha <- par.[f]
        D. <- matrix(par.[f1.1:f12],p-1,k)
        D. <- apply(D.,2, function(j) c(-sum(j), j))
        D. <- apply(D.,2, function(j) j-sum(j)/p)
        for (i in 1:k) {
            rss <- colSums((1/exp(alpha+D.[,i]))*(x-mu[,i])^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
        }
    },

    "VVI" = {
        alpha <- par.[f:f2]
        D. <- matrix(par.[f2.1:f22],p-1,k)
        D. <- apply(D.,2, function(j) c(-sum(j), j))
        D. <- apply(D.,2, function(j) j-sum(j)/p)
        for (i in 1:k) {
            rss <- colSums((1/exp(alpha[i]+D.[,i]))*(x-mu[,i])^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
        }
    },

    # here start the non-diagonal cases. main difference is the use
    # of backsolve() to calculate x^t Sigma^-1 x, works as follows:
    # assume Sigma = L D L^t, then Sigma^-1 = (L^t)^-1 D^-1 L^-1
    # y = L^-1 x  => x^t Sigma^-1 x = y^t D^-1 y
    # y = backsolve(L., x)

    "EEE" = {
        alpha <- par.[f]
        D. <- par.[f1.1:f11]
        D. <- c(-sum(D.), D.)
        D. <- D.-sum(D./p)
        invD <- (1/exp(alpha+D.))
        L. <- diag(1,p)
        L.[lower.tri(L., diag=FALSE)] <- par.[f11.1:f111]
        for (i in 1:k) {
            rss <- colSums(invD*backsolve(L.,(x-mu[,i]), upper.tri=FALSE)^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
        }
    },

    "VEE" = {
        alpha <- par.[f:f2]
        D. <- par.[f2.1:f21]
        D. <- c(-sum(D.), D.)
        D. <- D.-sum(D./p)
        L. <- diag(1,p)
        L.[lower.tri(L., diag=FALSE)] <- par.[f21.1:f211]
        for (i in 1:k) {
            rss <- colSums((1/exp(alpha[i]+D.))*backsolve(L., (x-mu[,i]), upper.tri=FALSE)^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
        }
    },

    "EVV" = {
        alpha <- par.[f]
        D. <- matrix(par.[f1.1:f12],p-1,k)
        D. <- apply(D.,2, function(j) c(-sum(j), j))
	    D. <- apply(D.,2, function(j) j-sum(j)/p)
        L.temp <- matrix(par.[f12.1:f121],p*(p-1)/2,k)
        for (i in 1:k) {
            L. <- diag(1,p)
            L.[lower.tri(L., diag=FALSE)] <- L.temp[,i]
            rss <- colSums((1/exp(alpha+D.[,i]))*backsolve(L., (x-mu[,i]), upper.tri=FALSE)^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha+l2pi)+rss))
        }
    },

    "VVV" = {
        alpha <- par.[f:f2]
        D. <- matrix(par.[f2.1:f22],p-1,k)
        D. <- apply(D.,2, function(j) c(-sum(j), j))
        D. <- apply(D.,2, function(j) j-sum(j)/p)
        invalpha <- (1/exp(rep(alpha, each=p)+D.))
        L.temp <- matrix(par.[f22.1:f221],p*(p-1)/2,k)
        L. <- diag(1,p)
        for (i in 1:k) {
            L.[lower.tri(L., diag=FALSE)] <- L.temp[,i]
            rss <- colSums(invalpha[,i]*backsolve(L., (x-mu[,i]), upper.tri=FALSE)^2)
            invl <- invl+w[i]*exp(-0.5*(p*(alpha[i]+l2pi)+rss))
        }
    },
    ## otherwise
    stop("invalid model:", model))

    ## return  sum_{i=1}^n log( f(x_i) ) :
    sum(log(invl))
}


#' log-likelihood function relying on mvtnorm function
#'
#' \code{llmvtnorm} returns scalar value of log-likelihood
#'
#' @param par. parameter vector as calculated by nMm2par
#' @param x matrix of samples
#' @param k number of cluster
#' @param trafo transformation of weights
#' @param model assumed model of the distribution
#'
#' @export
llmvtnorm <- function(par., x, k,
              trafo=c("clr1","logit"),
              model=c("EII","VII","EEI","VEI","EEE",
                  "VEE","EVI","VVI","EVV","VVV")
              )
{
    stopifnot(is.matrix(x),
              length(k <- as.integer(k)) == 1, k >= 1)
    trafo <- match.arg(trafo)
    model <- match.arg(model)
    p <- ncol(x)

    nmm <- par2nMm(par., p, k, trafo=trafo, model=model)
    ## FIXME (speed!):  dmvnorm(*, sigma= S) will do a chol(S) for each component
    ## -----  *instead* we already have LDL' and  chol(S) = sqrt(D) L' !!
    ## another par2*() function should give L and D, or from that chol(Sagma), rather than Sigma !
    w <- nmm$w
    mu <- nmm$mu
    sig <- nmm$Sigma
    y <- 0
    for (i in 1:k) {
        y <- y + w[i]*mvtnorm::dmvnorm(x,mean=mu[,i],sigma=sig[,,i])
    }
    sum(log(y))
}
