#### fit function for normal mixture samples


manyMLE <- function(x, k, models=1:10, 
                   trafo=c("clr1", "logit"),
                   ll = c("nmm", "mvt"),
                   savdir=NULL, name=NULL,
                ... ) {
    k <- as.integer(k)
    stopifnot(is.numeric(x),
              is.vector(models), length(models) <= 10,
              0 < models, models <= 10,
              is.integer(k), 0 < k)
    if (!is.null(savdir)) {stopifnot(dir.exists(savdir))}
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

    ret <- list(nMm=norMmixval, nMmtime=norMmixtime, k=k, 
                models=m, n=n, p=p, x=x)
    class(ret) <- "manyMLE"

    if (!is.null(savdir)) {
        savename <- ifelse(is.null(name), 
                           sprintf("fit_dim=%d_n=%d.rds", p, n),
                           paste0(name, ".rds"))
        saveRDS(list(fit=ret), file=file.path(normalizePath(savdir),savename))
    }
    ret
}


nobs.manyMLE <- function(object, ...) {object$n}


logLik.manyMLE <- function(object, ...) {
    ## returns log-likelihood of manyMLE object
    stopifnot(inherits(object, "manyMLE"))
    k <- object$k
    models <- object$models

    val <- matrix(0, length(k), length(models))
    rownames(val) <- k
    colnames(val) <- models

    for (i in seq_along(k)) {
        for (j in seq_along(models)) {
            nm <- object$nMm[i,j][[1]]
            # need to catch errors, if nm is string return NA
            val[i,j] <- ifelse(is.character(nm[[1]])&&length(nm)==2, 
                               NA, 
                               -nm$optr$value)
        }
    }

    #attributes(val) <- list(df=c(5,4,3), nobs=nobs(object))
    #class(val) <- c("logLik", "matrix")
    val
}

## DELETE_ME?
##displayError <- function(obj) {
##    stopifnot(inherits(obj, "manyMLE"))
##    k <- obj$k
##    models <- obj$models
##
##    for (i in seq_along(k)) {
##        for (j in seq_along(models)) {
##            nm <- obj$nMm[i,j][[1]]
##            if (is.character(nm[[1]])&&length(nm)==2) {
##                cat(k[i],models[j], "\t", paste(nm), "\n\n")
##            }
##        }
##    }
##}


npar.manyMLE <- function(object, ...) {
    k <- object$k
    p <- object$p
    models <- object$models

    val <- matrix(0, length(k), length(models))
    rownames(val) <- k
    colnames(val) <- models

    for (i in seq_along(k)) {
        for (j in seq_along(models)) {
            val[i,j] <- dfnMm(k[i],p,models[j])
        }
    }

    val
}


BIC.manyMLE <- function(object, ...) {
    n <- object$n
    k <- object$k
    models <- object$models
    npar <- npar(object)
    ll <- logLik(object)
    val <- npar*log(n) - 2*ll
    mi <- which.min(val)
    bestnMm <- object$nMm[mi][[1]]
    mirow <- mi%%length(k)
    micol <- ifelse(mirow>0, (mi%/%length(k))+1, mi%/%length(k))
    if (mirow==0) mirow <- length(k)
    mindex <- c(k[mirow],models[micol])
    list(val, best=mindex, bestnMm=bestnMm)
}


AIC.manyMLE <- function(object, ..., k = 2) {
    models <- object$models
    npar <- npar(object)
    ll <- logLik(object)
    val <- npar*k - 2*ll
    mi <- which.min(val)
    k <- object$k # overwriting the AIC k (typically = 2)
    mirow <- mi%%length(k)
    micol <- ifelse(mirow>0, (mi%/%length(k))+1, mi%/%length(k))
    if (mirow==0) mirow <- length(k)
    mindex <- c(k[mirow], models[micol])
    list(val, best=mindex)
}

## DELETE_ME?:
##cond <- function(obj) {
##    stopifnot(inherits(obj, "manyMLE")) # returns matrix of nobs/npar
##    k <- obj$k
##    models <- obj$models
##    val <- matrix(0, length(k), length(models))
##    rownames(val) <- k
##    colnames(val) <- models
##
##    for (i in seq_along(k)) {
##        for (j in seq_along(models)) {
##            nm <- obj$nMm[i,j][[1]]
##            # need to catch errors, if nm is string return NA
##            val[i,j] <- ifelse(is.character(nm[[1]])&&length(nm)==2, 
##                               NA, 
##                               nm$cond)
##        }
##    }
##
##    val
##}

## Needed in print() method:
extracttimes <- function(object, ...) {
    stopifnot(inherits(object, "manyMLE"))
    ti <- unlist(object$nMmtime)
    na <- names(ti)[1:5]
    co <- object$k
    mo <- object$models
    ti <- c(matrix(ti, ncol=5, byrow=TRUE))
    dn <- list(k=co, models=mo, proc_time=na)
    ## return
    structure(array(ti, lengths(dn), dimnames=dn),
              n = object$n, p = object$p, class = "manyMLE_time")
}
    

print.manyMLE <- function(x, ...) {
    n <- nobs(x)
    co <- x$k
    mo <- x$models
    dim <- x$p
    ti <- apply(extracttimes(x), 3, sum)
    bics <- BIC(x)

    cat("\nfitted normal mixture:\n",
        "dimension of dataset: \tvariables", dim, "\tobservations:", n, "\n",
        "fitted components and models: \n", co, "\n", mo, "\n")

    cat("total time: \t",ti, "\n")
    cat("\nbest fit:\t", bics[2][[1]], "\n",
        "logLik: \t", bics$bestnMm$optr$value, "\n")
    invisible(x)
}
