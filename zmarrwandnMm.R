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
		  weight = c(.2, .2, .6),
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
