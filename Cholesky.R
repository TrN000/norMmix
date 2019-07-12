### LDL Cholesky decomposition, not particularly efficient, but works for now


ldl <- function(mat){

	stopifnot( is.matrix(mat)&&is.numeric(mat)&&nrow(mat)==ncol(mat) )

	n <- ncol(mat)
	Diag <- rep(0,n)
	L <- matrix(0,n,n)
	mat.copy <- mat

	for (i in 1:n){
		Diag[i] <- mat[i,i]
		L[,i] <- mat[,i]/Diag[i]
		mat <- mat - Diag[i]*tcrossprod(L[,i],L[,i])
		mat[,i] <- rep(0,n)
		mat[i,] <- rep(0,n)
	}

	list(L=L, Diag=Diag)
}


# FF <- cbind( c(4.2,12,-16), c(12,37.07,-43), c(-16,-43,98.43) )
# FF
# ldl(FF)
# $L
#      [,1] [,2] [,3]
# [1,]    1    0    0
# [2,]    3    1    0
# [3,]   -4    5    1
# 
# $Diag
# [1] 4 1 9
# 

# system.time(  {for (i in 1:10000){ldl(FF)}}  )
#    user  system elapsed 
#   0.244   0.000   0.246 

# with isSymmetric
# system.time(  {for (i in 1:10000){ldl(FF)}}  )
#    user  system elapsed 
#   2.236   0.000   2.235 


# system.time(  {for (i in 1:1000000){ldl(FF)}}  )
#    user  system elapsed 
#  24.404   0.008  24.434 


