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


    bic <- matrix(0, length(k), length(m))
    colnames(bic) <- m
    rownames(bic) <- k

    aic <- matrix(0, length(k), length(m))
    colnames(aic) <- m
    rownames(aic) <- k

    for (j in 1:length(k)) {
        for (i in m) {
            nMm <- norMmixMLE(x,k[j],trafo=tr,model=i,...)
            bic[j,i] <- nMm$parlen*log(n) - 2*nMm$optr$value
            aic[j,i] <- nMm$parlen*2 - 2*nMm$optr$value
        }
    }

    indexb <- bic %in% sort(bic, decreasing=FALSE)[1:3]
    topbic <- which(matrix(indexb, dim(bic)[1], dim(bic)[2]), arr.ind=T)
    topbic <- cbind(topbic, bic[topbic])
    topbic[,1] <- k[topbic[,1]]
    topbic[,2] <- m[topbic[,2]]

    indexa <- aic %in% sort(aic, decreasing=FALSE)[1:3]
    topaic <- which(matrix(indexa, dim(aic)[1], dim(aic)[2]), arr.ind=T)
    topaic <- cbind(topaic, aic[topaic])
    topaic[,1] <- k[topaic[,1]]
    topaic[,2] <- m[topaic[,2]]

    list(BIC=bic, AIC=aic, topbic=topbic, topaic=topaic)

}

BIC.norMmixfit <- function(obj) {
    stopifnot(inherits(obj, "norMmixfit"))
    
    parlen <- obj$parlen
    n <- obj$n

    bic <- parlen*log(n) - 2*obj$optr$value
}

AIC.norMmixfit <- function(obj) {
    stopifnot(inherits(obj, "norMmixfit"))
    
    parlen <- obj$parlen

    bic <- parlen*2 - 2*obj$optr$value
}
	
    
