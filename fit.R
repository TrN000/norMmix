### fit function for normal mixture samples

#' @include norMmixMLE.R


fit.norMmix <- function(x, k, models=1:10, 
                        trafo=c("clr1", "logit"),
                        ll = c("nmm", "mvt"),
                    ...
                        )
{
    k <- as.integer(k)
    stopifnot(is.numeric(x),
              is.vector(models), length(models) <= 10,
              0 < models, models <= 10,
              is.integer(k), 0 < k)
    n <- nrow(x)
    p <- ncol(x)
    ll <- match.arg(ll)
    trafo <- match.arg(trafo)
    m <- c("EII","VII","EEI","VEI","EVI",
           "VVI","EEE","VEE","EVV","VVV")
    m <- m[models]

    norMmixval <- vector("list", length(m) * length(k))
    norMmixtime <- vector("list", length(m) * length(k))
    dim(norMmixval) <- c(length(k), length(m))
    dim(norMmixtime) <- c(length(k), length(m))

    for (j in seq_along(k)) {
        for (i in seq_along(m)) {
            st <- system.time(
                nMm <- tryCatch(norMmixMLE(x, k[j], model=m[i],
                                           ll=ll, trafo=trafo, ...), 
                                error = identity)
                )
            norMmixval[[j,i]] <- nMm
            norMmixtime[[j,i]] <- st
        }
    }

    rownames(norMmixval) <- k
    colnames(norMmixval) <- m
    rownames(norMmixtime) <- k
    colnames(norMmixtime) <- m

    ret <- list(nMm=norMmixval, nMmtime=norMmixtime, k=k, models=m, n=n, p=p, x=x)
    class(ret) <- c("fittednorMmix", "norMmix")
    ret
}

logLik.fittednorMmix <- function(object, ...)
{
    ## returns log-likelihood of fittednorMmix object
    stopifnot(inherits(object, "fittednorMmix"))
    k <- object$k
    models <- object$models

    val <- matrix(0, length(k), length(models))
    rownames(val) <- k
    colnames(val) <- models

    for (i in seq_along(k)) {
        for (j in seq_along(models)) {
            nm <- object$nMm[i,j][[1]]
            # need to catch errors, if nm is string return NA
            val[i,j] <- ifelse(is.character(nm[[1]])&&length(nm)==2, NA, -nm$optr$value)
        }
    }

    val
}


displayError.fittednorMmix <- function(obj)
{
    stopifnot(inherits(obj, "fittednorMmix"))
    k <- obj$k
    models <- obj$models

    for (i in seq_along(k)) {
        for (j in seq_along(models)) {
            nm <- obj$nMm[i,j][[1]]
            if (is.character(nm[[1]])&&length(nm)==2) cat(k[i],models[j], "\t", paste(nm), "\n\n")
        }
    }
}

parlen.fittednorMmix <- function(obj)
{
    stopifnot(inherits(obj, "fittednorMmix"))

    k <- obj$k
    p <- obj$p
    models <- obj$models

    val <- matrix(0, length(k), length(models))
    rownames(val) <- k
    colnames(val) <- models

    for (i in seq_along(k)) {
        for (j in seq_along(models)) {
            val[i,j] <- parlen(k[i],p,models[j])
        }
    }

    val
}


BIC.fittednorMmix <- function(object, ...)
{
    stopifnot(inherits(object, "fittednorMmix"))

    n <- object$n
    k <- object$k
    models <- object$models
    parlen <- parlen.fittednorMmix(object)
    ll <- logLik.fittednorMmix(object)
    val <- parlen*log(n) - 2*ll
    mi <- which.min(val)
    bestnMm <- object$nMm[mi][[1]]
    mirow <- mi%%length(k)
    micol <- ifelse(mirow>0, (mi%/%length(k))+1, mi%/%length(k))
    if (mirow==0) mirow <- length(k)
    mindex <- c(k[mirow],models[micol])

    list(val, best=mindex, bestnMm=bestnMm)
}

AIC.fittednorMmix <- function(object, ..., k = 2)
{
    stopifnot(inherits(object, "fittednorMmix"))

    ## n <- object$n
    models <- object$models
    parlen <- parlen.fittednorMmix(object)
    ll <- logLik.fittednorMmix(object)
    val <- parlen*k - 2*ll
    mi <- which.min(val)
    k <- object$k # overwriting the AIC k (typically = 2)
    mirow <- mi%%length(k)
    micol <- ifelse(mirow>0, (mi%/%length(k))+1, mi%/%length(k))
    if (mirow==0) mirow <- length(k)
    mindex <- c(k[mirow], models[micol])
    list(val, best=mindex)
}

cond.fittednorMmix <- function(obj)
{
    ## returns log-likelihood of fittednorMmix object

    stopifnot(inherits(obj, "fittednorMmix"))

    k <- obj$k
    models <- obj$models
    val <- matrix(0, length(k), length(models))
    rownames(val) <- k
    colnames(val) <- models

    for (i in seq_along(k)) {
        for (j in seq_along(models)) {
            nm <- obj$nMm[i,j][[1]]
            # need to catch errors, if nm is string return NA
            val[i,j] <- ifelse(is.character(nm[[1]])&&length(nm)==2, NA, nm$cond)
        }
    }

    val
}

