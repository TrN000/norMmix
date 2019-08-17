### fit function for normal mixture samples

#' @include norMmixMLE.R



fit.norMmix <- function(x, k=1:10, models=1:10, trafo=c("clr1","logit"),...) {
    stopifnot(is.numeric(x))
    stopifnot(is.vector(models), length(models)<=10, 
          isTRUE(all(models<=10)), isTRUE(all(models>0)))
    n <- nrow(x)
    p <- ncol(x)

    tr <- match.arg(trafo)

    m <- c("EII","VII","EEI","VEI","EVI",
        "VVI","EEE","VEE","EVV","VVV")
    m <- m[models]

    norMmixval <- list()


    for (j in 1:length(k)) {
        for (i in m) {
            nMm <- tryCatch(norMmixMLE(x,k[j],trafo=tr,model=i,...), error = function(e) paste("error",eval.parent(i,n=2),eval.parent(j,n=2)))
            norMmixval[[paste0(j,i)]] <- nMm
        }
    }

    ret <- list(nMm=norMmixval, k=k, models=m)
    class(ret) <- c("fittednorMmix", "norMmix")
    ret

}

BIC.fittednorMmix <- function(obj) {
    stopifnot(inherits(obj, "fittednorMmix"))

    k <- obj$k
    models <- obj$models

    


}

AIC.fittednorMmix <- function(obj) {
    stopifnot(inherits(obj, "fittednorMmix"))
    
}
	
    
