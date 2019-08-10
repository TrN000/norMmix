#### MLE for norMmix objects

#' @include norMmixEM.R





{}


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
               ini=c("cla","mcl"),
               maxiter=100, trace=2, tol=sqrt(.Machine$double.eps),
               samples=10,
	       ...
               ) {
    
    # 1. san check call
    # 2. prep nMm obj
    # 3. apply optim
    # 4. return


    # 1.
    trafo <- match.arg(trafo)
    model <- match.arg(model)

    stopifnot(is.numeric(x), is.numeric(k), (n <- nrow(x))>1)

    p <- ncol(x)
    k <- as.integer(k)
    

    ## init tau

    tau <- switch(ini,

    #init tau using clara
        "cla" = {
            clus <- cluster::clara(x=x, k, rngR=T, pamLike=T, samples=samples)
            index <- clus$clustering
            tau <- matrix(0,n,k)
            tau[cbind(1:n,index)] <- 1
            tau
            },

    # clustering using MBAhc from the mclust package
        "mcl" = {
            mclclus <- mclust::hcVVV(x)
            mclindex <- mclust::hclass(mclclus, k)
            mcltau <- matrix(0,n,k)
            mcltau[cbind(1:n,mclindex)] <- 1
            mcltau
            }
        )



    # 2.

    # one m-step
    nMm.temp <- mstep.nMm(x, tau,k)
    # create par. vector out of m-step
    initpar. <- nMm2par(obj=nMm.temp, trafo=trafo, model=model)

    #degrees of freedom
    parlen <- length(initpar.)
    

    # 3.

    # define function to optimize as negative log-lik
    # also reduces the number of arguments to par.
    neglogl <- function(par.) {
        -llmvtnorm(par.,x=x,p=p,k=k,trafo=trafo,model=model)
        }

    control <- list(maxit=maxiter, reltol = tol,
                    trace=(trace > 0), REPORT= pmax(1, 10 %/% trace),
    		    ...)

    optr <- optim(initpar., neglogl, method = "BFGS", control=control)

    optr$value <- -optr$value
    ### should be changed


    # 4.

    nMm <- par2nMm(optr$par, p, k, trafo=trafo, model=model)

    ret <- list(norMmix=nMm, optr=optr, parlen=parlen, n=n, mstep=nMm.temp, ini=initpar.)


    class(ret) <- "norMmixfit"

    return(ret)
}
