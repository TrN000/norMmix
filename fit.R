### fit function for normal mixture samples

#' @include norMmixMLE.R



fit.norMmix <- function(x, k=1:10, models=1:10, 
                        trafo=c("clr1","logit"),ll = c("nmm", "mvt"),
                    ...
                        )
{
    stopifnot(is.numeric(x),
              is.vector(models), length(models) <= 10,
              0 < models, models <= 10)
    n <- nrow(x)
    p <- ncol(x)

    trafo <- match.arg(trafo)
    ll <- match.arg(ll)

    m <- c("EII","VII","EEI","VEI","EVI",
           "VVI","EEE","VEE","EVV","VVV")
    m <- m[models]

    norMmixval <- list()


    for (j in 1:length(k)) {
        for (i in m) {
            nMm <- tryCatch(nMm <- norMmixMLE(x,k[j],trafo=trafo,model=i,ll=ll,...), error = identity)
            norMmixval[[paste0(i,j)]] <- nMm
        }
    }

    dim(norMmixval) <- c(length(m), length(k))
    norMmixval <- t(norMmixval)
    rownames(norMmixval) <- k
    colnames(norMmixval) <- m


    ret <- list(nMm=norMmixval, k=k, models=m, n=n, p=p, x=x)
    class(ret) <- c("fittednorMmix", "norMmix")
    ret

}

logLik.fittednorMmix <- function(obj)
{
    ## returns log-likelihood of fittednorMmix object

    stopifnot(inherits(obj, "fittednorMmix"))

    k <- obj$k
    models <- obj$models

    val <- matrix(0, length(k), length(models))
    rownames(val) <- k
    colnames(val) <- models

    for (i in k) {
        for (j in models) {
            nm <- obj$nMm[i,j][[1]]
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

    for (i in k) {
        for (j in models) {
            nm <- obj$nMm[i,j][[1]]
            if (is.character(nm[[1]])&&length(nm)==2) cat(i,j, "\t", paste(nm), "\n\n")
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

    for (i in k) {
        for (j in models) {
            val[i,j] <- parlen(i,p,j)
        }
    }

    val
}


BIC.fittednorMmix <- function(obj) 
{
    stopifnot(inherits(obj, "fittednorMmix"))

    n <- obj$n
    k <- obj$k
    models <- obj$models
    parlen <- parlen.fittednorMmix(obj)
    ll <- logLik.fittednorMmix(obj)
    val <- parlen*log(n) - 2*ll
    mi <- which.min(val)
    bestnMm <- obj$nMm[mi][[1]]
    mirow <- mi%%length(k)
    micol <- ifelse(mirow>0, (mi%/%length(k))+1, mi%/%length(k))
    if (mirow==0) mirow <- length(k)
    mindex <- c(k[mirow],models[micol])
    
    list(val, best=mindex, bestnMm=bestnMm)
}

AIC.fittednorMmix <- function(obj) 
{
    stopifnot(inherits(obj, "fittednorMmix"))

    n <- obj$n
    k <- obj$k
    models <- obj$models
    parlen <- parlen.fittednorMmix(obj)
    ll <- logLik.fittednorMmix(obj)
    val <- parlen*2 - 2*ll
    mi <- which.min(val)
    mirow <- mi%%length(k)
    micol <- ifelse(mirow>0, (mi%/%length(k))+1, mi%/%length(k))
    if (mirow==0) mirow <- length(k)
    mindex <- c(k[mirow],models[micol])
    
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

    for (i in k) {
        for (j in models) {
            nm <- obj$nMm[i,j][[1]]
            # need to catch errors, if nm is string return NA
            val[i,j] <- nm$cond
        }
    }

    val
}

