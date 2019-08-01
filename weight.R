
clr1 <- function(w) {

    stopifnot(is.numeric(w) ,all(w >= 0), all.equal(sum(w),1))

    # calculate clr1
    ln <- log(w)
    # return:
    ln[-1L] - mean(ln)
}


clr1inv <- function(p) {

    if (length(p)==0) {return(c(1))}
    stopifnot(is.numeric(p))
    # calc weights
    p1 <- exp(c(-sum(p), p))
    p1/sum(p1)
}


logit <- function(e) {

    stopifnot(is.numeric(e) ,all(e >= 0), all.equal(sum(e),1))

    qlogis(e[-1L])
}

logitinv <- function(e) {

    if (length(e)==0) {return(c(1))}
    stopifnot(is.numeric(e))

    e<- plogis(e)
    if ((sp. <- sum(e)) > 1) {
        stop("weights add to greater than 1")
    }
    
    w <- c((1-sp.), e)
}

