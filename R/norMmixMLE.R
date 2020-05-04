#### MLE for norMmix objects


## Compute clara()'s sampsize [ 'ss' := sampsize ;  'L' := Log ]
## to be used as argument e.g., of norMmixMLE()
ssClara2kL <- function(n,k, p)
    pmin(n,
         pmax(40, round(10*log(n))) + round(2*k*pmax(1, log(n*p))))

claraInit <- function(x, k,
                      samples = 128,
                      sampsize = ssClara2kL, # a number *or* function(n,k,p)
                      trace = 0)
{
    n <- nrow(x)
    p <- ncol(x)
    if(is.function(sampsize)) sampsize <- sampsize(n,k,p)
    stopifnot(length(sampsize) == 1L, sampsize >= 1)
    clus <- clara(x, k, rngR=TRUE, pamLike=TRUE, medoids.x=FALSE,
                  samples=samples, sampsize=sampsize, trace=trace)
    ## return index
    clus$clustering
}

## clustering using hc() from the mclust package
mclVVVinit <- function(x, k, ...) {
    mclclus <- hcVVV(x, ...)
    ## return index
    hclass(mclclus, k)
}

# Maximum likelihood Estimation for normal mixture models
#
# norMmixMLE returns fitted nMm obj
#
# Uses clara() and one M-step from EM-algorithm to initialize parameters
# after that uses general optimizer optim() to calculate ML.
#
# x:     sample matrix
# k:     number of components
# model: model to be assumed
norMmixMLE <- function(
               x, k,
               model = c("EII","VII","EEI","VEI","EVI",
                         "VVI","EEE","VEE","EVV","VVV"),
               initFUN = claraInit,
               ll = c("nmm", "mvt"),
               keep.optr = TRUE, keep.data = keep.optr,
               ## epsilon = 1e-10,
               method = "BFGS", maxit = 100, trace = 2,
               optREPORT=10, reltol = sqrt(.Machine$double.eps),
	       ... )
{
    # 1. san check call
    # 2. prep nMm obj
    # 3. apply optim
    # 4. return

    # 1.
    model <- match.arg(model)
    ll <- match.arg(ll)

    if(!is.matrix(x)) x <- data.matrix(x) # e.g. for data frame
    stopifnot(is.numeric(x), length(k <- as.integer(k)) == 1,
              k >= 1, (n <- nrow(x)) >= k,
              is.function(initFUN))
    p <- ncol(x)

    index <- initFUN(x,k)

    tau <- matrix(0, n,k)
    tau[cbind(1:n, index)] <- 1

    # 2.

    mcl.mstep <- switch(model,
        "EII" = mclust::mstepEII(x, tau),
        "VII" = mclust::mstepVII(x, tau),
        "EEI" = mclust::mstepEEI(x, tau),
        "VEI" = mclust::mstepVEI(x, tau),
        "EVI" = mclust::mstepEVI(x, tau),
        "VVI" = mclust::mstepVVI(x, tau),
        "EEE" = mclust::mstepEEE(x, tau),
        "VEE" = mclust::mstepVEE(x, tau),
        "EVV" = mclust::mstepEVV(x, tau),
        "VVV" = mclust::mstepVVV(x, tau),

        stop("error in mstep in norMmixMLE; model =", model)
    )

    stopifnot(is.list( par <- mcl.mstep$parameters ))
    nMm.temp <- norMmix(par$mean, Sigma = par$variance$sigma, weight = par$pro,
                        model = mcl.mstep$modelName)

    # create par. vector out of m-step
        #nMm.temp <- forcePositive(nMm.temp, eps0=epsilon)
    initpar. <- nMm2par(obj=nMm.temp, model=model,  meanFUN=mean)
    # save degrees of freedom
    npar <- length(initpar.)

    # 3.

    if(ll == "nmm") tx <- t(x)
    # define function to optimize as negative log-lik
    # also reduces the number of arguments to par.
    neglogl <-
        switch(ll,
               "nmm" = function(P) -llnorMmix(P, tx=tx, k=k, model=model),
               ## max(-10^300, -llnorMmix) for both
               "mvt" = function(P) -llmvtnorm(P,  x= x, k=k, model=model),
               stop("error selecting neglogl")) # impossible as have ll <- match.arg(.)

    control <- list(maxit=maxit, reltol=reltol,
                    trace = (trace > 0), REPORT= optREPORT,
    		    ...)
    optr <- optim(initpar., neglogl, method=method, control=control)

    ## 4.  return ()

    # newer data structure, but still not inheriting from "norMmix" but rather with 'x$norMmix':
    structure(class = "norMmixMLE",
              list(norMmix = par2nMm(optr$par, p, k, model=model)
                 , npar = npar
                 , nobs = n
                 , cond = parcond(x, k=k, model=model)
                 , optr = if(keep.optr) optr
                 , x    = if(keep.data) x))

    ##new data structure
    ##ret <- structure(par2nMm(optr$par, p, k, model=model),
    ##                 class=c("norMmixMLE","norMmix"),
    ##                 nobs=n,
    ##                 npar=npar,
    ##                 cond=parcond(x, k=k, model=model))
    ##if (keep.optr) attr(ret, "optr") <- optr
    ##if (keep.data) attr(ret, "x") <- x
    ## ret
}

## superfluous:
## is.norMmixMLE <- function(x) inherits(x, "norMmixMLE")


logLik.norMmixMLE <- function(object, ...) {
    structure( - object$optr$value
            , class = "logLik"
            , df    = object$npar
            , nobs  = object$nobs)
}

nobs.norMmixMLE <- function(object, ...) object$nobs

npar.norMmixMLE <- function(object, ...) object$npar ## ==!?== npar(object$norMmix)

##' format 'x' with names into  "<name1> = <x1>, <name2> = <x2>, ..."
named2char <- function(x, sep = " = ", collapse = ", ")
    paste(names(x), x, sep=sep, collapse=collapse)

print.norMmixMLE <- function(x, ...) {
    cat("'norMmixMLE' normal mixture MLE fit;  fitted 'norMmix' normal mixture:\n")
    print.norMmix(x$norMmix, ...)
    cat("log-likelihood:", x$logLik, "\n\n",
        "nobs\t npar\t nobs/npar\n",
        x$nobs, "\t", x$npar, "\t", x$cond, "\n")
    if(!is.null(optr <- x$optr))
        cat("\n optim() counts:", named2char(optr$counts),"\n")
    invisible(x)
}
