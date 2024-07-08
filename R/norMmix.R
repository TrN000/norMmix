## norMmix data structure


## Auxiliary function evals to TRUE if x is sym pos def array, otherwise character string with msg
okSigma <- function(Sig, tol1 = 1000 * .Machine$double.eps, tol2 = 1e-10) {
    if (!is.numeric(Sig) || !is.array(Sig) || length(d <- dim(Sig)) != 3) {
        "is not a numeric array of rank 3"
    } else if (prod(d) == 0) {
        "at least one of the dimensions is zero"
    } else { ## all d[.] >= 1
        k <- d[[3]]
        p <- d[[1]]
        if (d[[2]] != p) {
            "Sigma matrix dimension is not square"
        } else {
            if(p >= 2) for (i in 1:k) {
                Si <- Sig[, , i]
                if (!isSymmetric(Si, tol = tol1)) {
                    return(paste0("Sigma[ , ,", i, "] is not symmetric"))
                }
                ## else
                ev <- eigen(Si, only.values = TRUE)$values
                if (any(ev < -tol2)) {
                    return(paste0("Sigma[ , ,", i, "] is not positive semi-definite"))
                }
            }
            TRUE
        }
    }
}


# Constructor for nMm 'objects'

norMmix <- function(mu,
                    Sigma = NULL,
                    weight = rep(1 / k, k),
                    name = NULL,
                    model = c(
                        "EII", "VII", "EEI", "VEI", "EVI",
                        "VVI", "EEE", "VEE", "EVV", "VVV"
                    )) {
    ## Purpose: constructor for 'norMmix' (multivariate normix)
    ## --------------------------------------------------------
    ## Arguments:  ===> ../man/norMmix.Rd
    ##
    ##    mu: matrix with columns as vector means of dist.
    ##    Sigma:     option 0: default, generates all identity
    ##            covariance mat.
    ##        option 1: scalar, generates EII dist
    ##        option 2: vector of length k, generates VII
    ##            distribution
    ##        option 3: array of dimension p x p x k.
    ##            covariance matrices of distributions
    ##    weight: vector of length k, sums to 1
    ##    name: name attribute
    ##    model: type of distribution VVV, IVV etc.

    ### FIXME:  1) model is *NOT* checked
    ### -----   2) For simple 'model', simple 'Sigma = ..' make sense model-*specifically* !
    ## ==> unfinished  guessModel  below

    ## --------------------------------------------------------
    ## Value: returns objext of class 'norMmix'
    ## --------------------------------------------------------
    ## Author: nicolas trutmann, Date:2019-06-19

    stopifnot(is.numeric(mu))
    guessModel <- missing(model) # if not specified, may guess it from {p, k, Sigma}

    if (!is.matrix(mu)) mu <- as.matrix(mu) # p x 1  typically
    # TODO: if no Sigma is supplied might tacitly mistake px1 and 1xp; throw error
    p <- nrow(mu) # p = dimension
    k <- ncol(mu) # k = number of components

    dS <- c(p, p, k) # == dim(Sigma) {in full VVV  parametrization !}
    if (is.null(Sigma)) {
        Sigma <- 1
    } # (option 0)
    else {
        stopifnot(is.numeric(Sigma))
    }
    isArr <- is.array(Sigma) # if not array, is also not matrix
    if (!isArr && length(Sigma) == 1) {
        Sigma <- array(diag(Sigma, p), dS)
    } else if (!isArr && length(Sigma) == k) {
        Sigma <- array(vapply(Sigma, function(s) diag(s, p), diag(p)), dS)
    } else if (!(isArr && all(dim(Sigma) == c(p, p, k)))) {
        stop("'Sigma' not among recognized formats")
    }
    if (!isTRUE(m <- okSigma(Sigma))) stop(m)

    # inspect weight
    stopifnot(is.numeric(weight))
    if (length(weight) != k) {
        stop("weight is not of length k")
    }
    if (!(all(weight >= 0) && (abs(sum(weight) - 1) < 1000 * .Machine$double.eps))) {
        stop("weight doesn't sum to 1 or isn't positive")
    }

    model <- if (missing(model)) "VVV" else match.arg(model)
    if (is.null(name)) {
        name <- sprintf("model \"%s\", G = %s", model, k)
    }
    structure(
        name = name,
        class = "norMmix",
        list(
            model = model,
            mu = mu, Sigma = Sigma, weight = weight,
            k = k # == length(weight) == ncol(mu)
            , dim = p # == nrow(mu)
        )
    )
}

npar.norMmix <- function(object, ...) {
    mo <- object$model
    k <- object$k
    p <- object$dim
    dfnMm(k, p, mo)
}


is.norMmix <- function(x) {
    inherits(x, "norMmix")
}

nor1toMmix <- function(object) {
    mu <- matrix(object[, 1], nrow = 1, byrow = TRUE)
    k <- length(object[, 2])
    sigma <- array(object[, 2], c(1, 1, 2))
    weights <- object[, 3]
    name <- attr(object, "name")

    norMmix(mu, Sigma = sigma, weight = weights, name = name, model = "VII")
}

## -- *not* used currently!
#' Corrects numerical error in case D in LDL' decomp is near-zero:
#' takes tolerance, norMmix obj; returns norMmix obj
# forcePositive <- function(nMm, eps0=1e-10) {
#    stopifnot(is.norMmix(nMm))
#
#    sig <- nMm$Sigma
#    D. <- apply(sig,3, function(j) ldl(j)$D )
#    L. <- apply(sig,3, function(j) ldl(j)$L)
#    k <- ncol(D.)
#    p <- nrow(D.)
#
#    eps <- eps0 * apply(D., 2, function(dj) max(abs(dj)))
#
#    for (i in 1:k) { # for each component
#        Di <- pmax(eps[i], D.[,i])
#        ## Sig := L D L' :
#        sig[,,i] <- matrix(L.[,i],p,p) %*% diag(Di) %*% matrix(L.[,i],p,p,byrow=TRUE)
#    }
#
#    nMm$Sigma <- sig
#    nMm
# }


# rnorMmix
rnorMmix <- function(n, obj, index = FALSE, permute = TRUE) {
    ## Purpose: generates random values distributed by NMD
    ## -------------------------------------------------------------------
    ## Arguments:
    ##    n: number of p-dimensional observations desired
    ##  obj: of type norMmix
    ## -------------------------------------------------------------------
    ## Value: matrix  n x p (columns are vectors)
    ## -------------------------------------------------------------------
    ## Author: nicolas trutmann, Date:2019-06-21
    if (!inherits(obj, "norMmix")) stop("argument must be of class \"norMmix\"")

    mu <- obj$mu
    Sigma <- obj$Sigma
    weight <- obj$weight
    ## p <- obj$dim
    nj <- rmultinom(n = 1, size = n, prob = weight)
    ### FIXME: again this does  chol(Sigma_j)  j = 1...k .. when actually we could've rather stored  chol(Sigma) instead of Sigma
    a <- do.call(
        rbind,
        lapply(seq_along(nj), function(j) mvrnorm(n = nj[j], mu = mu[, j], Sigma = Sigma[, , j]))
    )
    if (index) {
        cl <- rep(seq_along(weight), times = nj)
        a <- cbind(cl, a)
    }
    if (permute) a[sample.int(n), ] else a
}


# density function for norMmix object
dnorMmix <- function(x, nMm) {
    stopifnot(
        is.norMmix(nMm), is.numeric(x),
        length(p <- nMm$dim) == 1, p >= 1,
        length(k <- nMm$ k) == 1, k >= 1
    )
    if (!is.matrix(x)) x <- as.matrix(x)
    stopifnot(ncol(x) == p)
    ## FIXME: Using dmvnorm() is slow as it needs to solve( Sigma_j ) for each  j = 1..k
    ret <- 0
    for (i in 1:k) {
        ret <- ret + nMm$weight[i] * mvtnorm::dmvnorm(x, mean = nMm$mu[, i], sigma = nMm$Sigma[, , i])
    }
    ret
}

dnorMmixL <- function(obj, x = NULL, log = FALSE, xlim = NULL, n = 511) {
    if (!inherits(obj, "norMmix")) stop("argument must be of class \"norMmix\"")
    if (is.null(x)) {
        if (is.null(xlim)) {
            xlim <- mean.norMmix(obj) + c(-3, 3) * sqrt(var.norMix(obj))
        }
        x <- seq(xlim[1], xlim[2], length = n)
    }
    list(x = x, y = dnorMmix(x, obj))
}


mean.norMmix <- function(x, ...) {
    if (!inherits(x, "norMmix")) stop("argument must be of class \"norMmix\"")
    drop(x$mu %*% x$weight)
}


var.norMix <- function(x, ...) {
    if (!inherits(x, "norMmix")) stop("argument must be of class \"norMmix\"")
    stop("NOT IMPLEMENTED YET") # TODO: figure out what variance is.
}


mods <- function() { ## just return
    c(
        "EII", "VII", "EEI", "VEI", "EVI",
        "VVI", "EEE", "VEE", "EVV", "VVV"
    )
}


print.norMmix <- function(x, ...) {
    name <- attr(x, "name")
    mo <- x$model
    we <- x$weight
    co <- x$k
    dim <- x$dim

    cat("norMmix object: \n")
    cat("multivariate normal mixture model with the following attributes:\n")
    cat(
        "name: \t\t", name, "\n",
        "model: \t\t", mo, "\n",
        "dimension:\t", dim, "\n",
        "components:\t", co, "\n"
    )
    cat(
        "weight of components",
        sort(signif(we, digits = 3), decreasing = TRUE), "\n"
    )
    invisible(x) # << standard for all "good citizen" print() methods
}
