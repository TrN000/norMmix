
clr1 <- function(w) {

	stopifnot(is.numeric(w) ,all(w >= 0), all.equal(sum(w),1))

	# calculate clr1
	ln <- log(w)
	res <- ln[-1L] - mean(ln)
	res
}


clr1inv <- function(p) {

	if (length(p)==0) {return(c(1))}
	stopifnot(is.numeric(p))
	# calc weights
	p1 <- exp(c(-sum(p), p))
	f <- sum(p1)
	pp <- p1/f
	pp
}
