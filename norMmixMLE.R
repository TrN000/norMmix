#### MLE for norMmix objects

#' @include norMmixEM.R

{}


##' Compute clara()'s sampsize [ 'ss' := sampsize ;  'L' := Log ]
##' to be used as argument e.g., of norMmixMLE()
##' @export
ssClaraL <- function(n,k, p) pmin(n, pmax(40, round(10*log(n))) + round(2*k*pmax(1, log(n*p))))


#' Maximum likelihood Estimation for normal mixture models
#'
#' \code{norMmixMLE} returns fitted nMm obj
#'
#' Uses clara() and one M-step from EM-algorithm to initialize parameters
#' after that uses general optimizer optim() to calculate ML.
#'
#' @param x sample matrix
#' @param k number of clusters
#' @param trafo transformation to be used
#' @param model model to be assumed
#'
#' @export
norMmixMLE <- function(
               x, k,
               trafo = c("clr1", "logit"),
               model = c("EII","VII","EEI","VEI","EVI",
                         "VVI","EEE","VEE","EVV","VVV"),
               ini = c("clara", "mclVVV"),
               ll = c("nmm", "mvt"),
               method = "BFGS", maxit = 100, trace = 2, reltol = sqrt(.Machine$double.eps),
               samples = 128,
               sampsize = ssClaraL,
               traceClara = 0,
	       ...
               )
{
    # 1. san check call
    # 2. prep nMm obj
    # 3. apply optim
    # 4. return


    # 1.
    trafo <- match.arg(trafo)
    model <- match.arg(model)
    ll <- match.arg(ll)

    if(!is.matrix(x)) x <- data.matrix(x) # e.g. for data frame
    stopifnot(is.numeric(x), length(k <- as.integer(k)) == 1, (n <- nrow(x)) > 1)
    p <- ncol(x)

    ## init tau : index <- <clustering>
    switch(ini,
           ## init tau using clara
        "clara" = {
            if(is.function(sampsize)) sampsize <- sampsize(n,k,p)
            stopifnot(length(sampsize) == 1L, sampsize >= 1)
            clus <- clara(x, k, rngR=TRUE, pamLike=TRUE, medoids.x=FALSE,
                          samples=samples, sampsize=sampsize, trace=traceClara)
            index <- clus$clustering
        },

        ## clustering using hc() from the mclust package
        "mclVVV" = {
            mclclus <- hcVVV(x)
            index <- hclass(mclclus, k)
        },
        stop("invalid 'ini':", ini))

    tau <- matrix(0, n,k)
    tau[cbind(1:n, index)] <- 1

    # 2.

    # one M-step  (TODO: mstep() could *depend* on 'model'; currently does "VVV")
    nMm.temp <- mstep.nMm(x, tau)
    # create par. vector out of m-step

    nMm.temp <- forcePositive(nMm.temp)

    initpar. <- nMm2par(obj=nMm.temp, trafo=trafo, model=model, meanFUN=mean)
    #degrees of freedom
    parlen <- length(initpar.)


    # 3.

    # define function to optimize as negative log-lik
    # also reduces the number of arguments to par.
    neglogl <- switch(ll,
        "nmm" = function(par.) { -llnorMmix(par.,x=t(x),k=k,trafo=trafo,model=model) },
        "mvt" = function(par.) { -llmvtnorm(par.,x=x,k=k,trafo=trafo,model=model) },
        stop("error selecting neglogl") )

    control <- list(maxit=maxit, reltol=reltol,
                    trace=(trace > 0), REPORT= pmax(1, 10 %/% trace),
    		    ...)
    optr <- optim(initpar., neglogl, method=method, control=control)


    # 4.

    nMm <- par2nMm(optr$par, p, k, trafo=trafo, model=model)
    cond <- parcond(x, k=k, model=model)

    ret <- list(norMmix=nMm, optr=optr, parlen=parlen, n=n, mstep=nMm.temp, ini=initpar., cond=cond)

    class(ret) <- "norMmixfit"
    return(ret)
}
