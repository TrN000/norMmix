### fit function for normal mixture samples

#' @include norMmixMLE.R



fit.norMmix <- function(x, k=1:10, models=1:10, trafo=c("clr1","logit"),...) {
    stopifnot(is.numeric(x),
              is.vector(models), length(models) <= 10,
              0 < models, models <= 10)
    n <- nrow(x)
    p <- ncol(x)

    trafo <- match.arg(trafo)

    m <- c("EII","VII","EEI","VEI","EVI",
           "VVI","EEE","VEE","EVV","VVV")
    m <- m[models]

    norMmixval <- list()


    for (j in 1:length(k)) {
        for (i in m) {
            nMm <- tryCatch(nMm <- norMmixMLE(x,k[j],trafo=trafo,model=i,...), error = function(e) paste("error",eval.parent(i,n=2),eval.parent(j,n=2)))
            norMmixval[[paste0(i,j)]] <- nMm
        }
    }

    ret <- list(nMm=norMmixval, k=k, models=m, n=n)
    class(ret) <- c("fittednorMmix", "norMmix")
    ret

}

logLik.fittednorMmix <- function(obj) {
    stopifnot(inherits(obj, "fittednorMmix"))

    k <- obj$k
    models <- obj$models


    val <- matrix(0, length(k), length(models))
    rownames(val) <- k
    colnames(val) <- models

    for (i in k) {
        for (j in models) {
            val[i,j] <- obj$nMm[[paste0(j,i)]]$optr$value
        }
    }

    val
}


BIC.fittednorMmix <- function(obj) {
    stopifnot(inherits(obj, "fittednorMmix"))

    n <- obj$n
    k <- obj$k
    models <- obj$models

    parlen <- matrix(0, length(k), length(models))
    rownames(parlen) <- k
    colnames(parlen) <- models

    for (i in k) {
        for (j in models) {
            parlen[i,j] <- obj$nMm[[paste0(j,i)]]$parlen
        }
    }

    ll <- logLik.fittednorMmix(obj)

    val <- parlen*log(n) - 2*ll

    val

}
