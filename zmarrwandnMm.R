### filename starts with 'z' so it gets loaded after norMmix.R

### After Steve Marron, adapted for multivariate normal mixtures


### 2D mixtures

#' @include norMmix.R


MW2nm1 <- norMmix(
		  name = "#1 gaussian",
		  mu = cbind( c(0,0)), 
		  Sigma = c(1),
		  model = "EII"
		  )

MW2nm1.2 <- norMmix(
		    name = "one component rotated",
		    mu = cbind( c(0,0) ),
		    Sigma = array(cbind(c(55,9), c(9,3)), c(2,2,1)),
		    model = "EVV"
		    )

MW2nm2 <- norMmix(
		  name = "#2 skewed",
		  mu = cbind( c(0,0), c(0.5,0), c(13/12,0)),
		  Sigma = c(1, (2/3), (5/9)),
		  weight = c(.6, .2, .2),
		  model = "VII"
		  )

MW2nm4 <- norMmix(
		  name = "#4 kurtotic",
		  mu = cbind( c(0,0), c(0,0)),
		  Sigma = c(1,.1),
		  weight = c(2/3, 1/3),
		  model = "VII"
		  )

MW2nm5 <- norMmix(
		  name = "#5 test5",
		  mu = cbind( c(0,0), c(0,0)),
		  Sigma = c(1,.1),
		  weight = c(2/3, 1/3),
		  model = "VII"
		  )

MW26 <- norMmix(
		name = "#6 test EEI",
		weight = c(0.5, 0.5),
		mu = cbind( c(0,0), c(3,3) ),
		Sigma = array(rep(diag(c(4,5)),2),c(2,2,2)),
		model = "EEI"
		)

MW27 <- norMmix(
		name = "#7 test VEI",
		weight = c(0.5, 0.5),
		mu = cbind( c(0,0), c(3,3) ),
		Sigma = array(c(diag(c(4,5)),diag(c(8,10))),c(2,2,2)),
		model = "VEI"
		)

MW28 <- norMmix(
		name = "#8 test EVI",
		weight = c(0.2, 0.2, 0.6),
		mu = cbind( c(0,0), c(1,1), c(-1,-1) ),
		Sigma = array(c(diag(c(2,9)),diag(c(9,2)),diag(c(3,6))),c(2,2,3)),
		model = "EVI"
		)

MW29 <- norMmix(
		name = "#9 test VVI",
		weight = c(0.5, 0.5),
		mu = cbind( c(0,0), c(3,3) ),
		Sigma = array(c(diag(c(4,5)),diag(c(7,11))),c(2,2,2)),
		model = "VVI"
		)

MW210 <- norMmix(
		name = "#10 test EEE",
		weight = c(0.5, 0.5),
		mu = cbind( c(0,0), c(3,3) ),
		Sigma = array(c( 1,3,3,11, 1,3,3,11 ),c(2,2,2)),
		model = "EEE"
		)

MW211 <- norMmix(
		name = "#11 test VEE",
		weight = c(0.5, 0.5),
		mu = cbind( c(0,0), c(3,3) ),
		Sigma = array(c( 1,3,3,11, 2,6,6,22 ),c(2,2,2)),
		model = "VEE"
		)

MW212 <- norMmix(
		name = "#12 test EVV",
		weight = c(0.5, 0.5),
		mu = cbind( c(0,0), c(3,3) ),
		Sigma = array(c( 1,3,3,11, 2,4,4,9 ),c(2,2,2)),
		model = "EVV"
		)

MW213 <- norMmix(
		name = "#13 test VVV",
		weight = c(0.5, 0.5),
		mu = cbind( c(0,0), c(30,30) ),
		Sigma = array(c( 1,3,3,11, 3,6,6,13 ),c(2,2,2)),
		model = "VVV"
		)



n2m <- function(obj) nMm2par(obj , trafo="clr1", obj$model)

par210 <- n2m(MW210)
par211 <- n2m(MW211)
par212 <- n2m(MW212)
par213 <- n2m(MW213)
par26 <- n2m(MW26)
par27 <- n2m(MW27)
par28 <- n2m(MW28)
par29 <- n2m(MW29)
par2nm1 <- n2m(MW2nm1)
par2nm1.2 <- n2m(MW2nm1.2)
par2nm2 <- n2m(MW2nm2)
par2nm4 <- n2m(MW2nm4)
par2nm5 <- n2m(MW2nm5)


ln2m <- function(obj) nMm2par(obj , trafo="logit", obj$model)


lpar210 <- ln2m(MW210)
lpar211 <- ln2m(MW211)
lpar212 <- ln2m(MW212)
lpar213 <- ln2m(MW213)
lpar26 <- ln2m(MW26)
lpar27 <- ln2m(MW27)
lpar28 <- ln2m(MW28)
lpar29 <- ln2m(MW29)
lpar2nm1 <- ln2m(MW2nm1)
lpar2nm1.2 <- ln2m(MW2nm1.2)
lpar2nm2 <- ln2m(MW2nm2)
lpar2nm4 <- ln2m(MW2nm4)
lpar2nm5 <- ln2m(MW2nm5)
