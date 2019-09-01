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
                "VVI","EEE","VEE","EVV","VVV"),
            meanFUN= mean
            ){

    #transferring values of obj to handier variables
    w <- obj$weight
    mu <- obj$mu
    sig <- obj$Sigma
    p <- obj$dim
    k <- obj$k

    trafo <- match.arg(trafo)
    model <- match.arg(model)

    av <- match.fun(meanFUN)

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
    stopifnot( isTRUE( all( apply(sig,3, function(j) (ldl(j)$D >= 0 )))))


    #output vector of parameter values

    c(
      w <- switch(trafo, #weights either logit or centered log ratio
        "logit" = logit(w),

            "clr1" = clr1(w),

        stop("Error in weight trafo, ",trafo)
            ),
      mu, #means
      Sigma <- switch(model, #model dependent covariance values
        "EII" = {
            D. <- apply(sig,3, function(j) ldl(j)$D)
            av(log(D.))
            },

        "VII" = {
            D. <- apply(sig,3, function(j) ldl(j)$D)
            apply(D.,2, function(j) av(log(j)))
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
            D. <- apply(D.,1, function(j) av(log(j)))
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
            D. <- apply(sig,3, function(j) ldl(j)$D)
            alpha <- apply(D.,2, function(j) av(log(j)))
            D. <- apply(D.,2, function(j) log(j)-av(log(j)))
            c(alpha, D.[-1,])
            },

        "EEE" = {
            D. <- apply(sig,3, function(j) log(ldl(j)$D))
            D. <- apply(D.,1, av)
            alpha <- av(D.)
            D. <- D.-av(D.)
            L. <- apply(sig,3, function(j) ld.(ldl(j)$L))
            L. <- matrix(L., p*(p-1)/2, k)
            L. <- apply(L.,1, function(j) av(j))
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

        stop("invalid argument in 'model'")
        )
    )
}


#' wrapper function for nMm objs in zmarrwandMm
#'
#' \code{n2p} returns same as nMm2par with clr1
#'
#' @export
#n2p <-
## these were in ./zmarrwandnMm.R :
#n2m <- # <- drop this name and rather use
nc2p <- function(obj) nMm2par(obj , trafo="clr1",  obj$model)

#ln2m <- # <- drop this, and rather use
nl2p <- function(obj) nMm2par(obj , trafo="logit", obj$model)



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
            trafo = c("clr1", "logit"),
            model = c("EII","VII","EEI","VEI","EEE",
                    "VEE","EVI","VVI","EVV","VVV")
            )
{
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

    f11 <- f1 + p -1 # end of D. if D. uniform and alpha uniform
    f12 <- f1 + p*k -k # end D. if D. var and alpha unif.
    f21 <- f2 + p -1 # end of D. if D. uniform and alpha variable
    f22 <- f2 + p*k -k # end of D. if D.var and alpha var

    f11.1 <- f11 +1L # start of L if alpha unif D unif
    f21.1 <- f21 +1L # start of L if alpha var D unif
    f12.1 <- f12 +1L # start of L if alpha unif D var
    f22.1 <- f22 +1L # start of L if alpha var D var

    f111 <- f11 + p*(p-1)/2 # end of L if alpha unif D unif
    f211 <- f21 + p*(p-1)/2 # end of L if alpha var D unif
    f121 <- f12 + k*p*(p-1)/2 # end of L if alpha unif D var
    f221 <- f22 + k*p*(p-1)/2 # end of L if alpha var D var

    #only important ones are f1.2, f1.3, f2.2, f2.3

    w.temp <- if(k==1) vector() else par.[1:(k-1)]
    w <- switch(trafo,
                "logit" = logitinv(w.temp),
                "clr1"  = clr1inv (w.temp),
                stop("invalid 'trafo'": trafo))

    mu <- matrix(par.[k:(k+p*k-1)], p, k)

### FIXME: Alternatively, instead of Sigma, compute  chol(Sigma) = D^{1/2} L'  as Sigma = LDL'

    Sigma <- switch(model,
    # diagonal cases
    "EII" = {
        lambda <- exp(par.[f])
        array( rep(diag(lambda, p),k), c(p,p,k) )
        },

    "VII" = {
        lambda <- exp(par.[f:f2])
        array(unlist(lapply( lambda, function(j) diag(j,p) )), c(p,p,k))
        },

    "EEI" = {
        lambda <- par.[f]
        D. <- par.[f1.1:f11]
        D. <- c(-sum(D.), D.)
        D. <- D.-mean(D.)
        array( rep(diag(exp(lambda+D.)),k), c(p,p,k) )
        },

    "VEI" = {
        lambda <- par.[f:f2]
        D. <- par.[f2.1:f21]
        D. <- c(-sum(D.), D.)
        D. <- matrix(D.+rep(lambda, each=p), p, k)
        D. <- exp(D.)
        array( apply(D.,2, diag), c(p,p,k))
        },

    "EVI" = {
        lambda <- par.[f]
        D. <- matrix(par.[f1.1:f12],p-1,k)
        D. <- apply(D., 2, function(j) c(-sum(j), j))
        D. <- exp(D.+lambda)
        array( apply(D.,2, diag), c(p,p,k))
        },

    "VVI" = {
        lambda <- par.[f:f2]
        D. <- matrix(par.[f2.1:f22],p-1,k)
        D. <- apply(D., 2, function(j) c(-sum(j), j))
        D. <- exp(D.+rep(lambda, each=p))
        array(apply(D.,2, diag), c(p,p,k))
        },

    # variable cases

    "EEE" = {
        lambda <- par.[f]
        D. <- par.[f1.1:f11]
        D. <- c(-sum(D.), D.)
        D. <- exp(D.+lambda)
        L. <- par.[f11.1:f111]
        A. <- dl.(D.,L.,p)
        sig <- array(rep(A., times=k), c(p,p,k))
        sig
        },

    "VEE" = {
        lambda <- par.[f:f2]
        D. <- par.[f2.1:f21]
        D. <- c(-sum(D.), D.)
        D. <- exp(matrix(D.+rep(lambda, each=p), p, k))
        f3 <- (p*(p-1)/2)
        L. <- par.[f21.1:f211]
        sig <- array(0, c(p,p,k))
        for (i in 1:k){
            sig[,,i] <- dl.(D.[,i],L.,p)
        }
        sig
        },

    "EVV" = {
        #par.[f:f12] <- exp(par.[f:f12])
        lambda <- par.[f]
        D. <- matrix(par.[f1.1:f12],p-1,k)
        D. <- apply(D., 2, function(j) c(-sum(j), j))
        D. <- exp(D.+lambda)
        f3 <- (p*(p-1)/2)
        L.temp <- matrix(par.[f12.1:f121],f3,k)
        sig <- array(0, c(p,p,k))
        for (i in 1:k) {
            sig[,,i] <- dl.(D.[,i],L.temp[,i],p)
        }
        sig
        },

    "VVV" = {
        #par.[f:f22] <- exp(par.[f:f22])
        lambda <- par.[f:f2]
        D. <- matrix(par.[f2.1:f22],p-1,k)
        D. <- apply(D., 2, function(j) c(-sum(j), j))
        D. <- exp(D.+rep(lambda,each=p))
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




parlen <- function(k,p,
                   model=c("EII","VII","EEI","VEI","EVI",
                           "VVI","EEE","VEE","EVV","VVV")
                   ) {

    stopifnot(is.numeric(k), is.numeric(p))
    model <- match.arg(model)

    w <- k-1
    mu <- p*k
    sig <- switch(model,
        
        "EII" = 1,

        "VII" = k,

        "EEI" = 1+ (p-1),

        "VEI" = k+ (p-1),

        "EVI" = 1+ k*(p-1),

        "VVI" = k+ k*(p-1),

        "EEE" = 1+ (p-1)+ p*(p-1)/2,

        "VEE" = k+ (p-1)+ p*(p-1)/2,

        "EVV" = 1+ k*(p-1)+ k*p*(p-1)/2,

        "VVV" = k+ k*(p-1)+ k*p*(p-1)/2
        )

    param <- w+mu+sig
    param
}


parcond <- function(x,
                    k,
                    model=c("EII","VII","EEI","VEI","EVI",
                            "VVI","EEE","VEE","EVV","VVV")
                    ) {

    n <- nrow(x)
    p <- ncol(x)
    model <- match.arg(model)
    pars <- parlen(k,p, model=model)

    n/pars
}
                    
