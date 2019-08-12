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


####----------------------------------------------------------------------------
## 3 dims
####----------------------------------------------------------------------------

MW31 <- norMmix(
    name = "#1 3d EII",
    weight = 1,
    mu = as.matrix(c(0,0,0)),
    Sigma = c(1),
    model = "EII"
    )

MW32 <- norMmix(
    name = "#2 3d VII",
    weight = c(0.2,0.2,0.2,0.2,0.2),
    mu = matrix(1:15, 3,5),
    Sigma = 1:5,
    model = "VII"
    )

MW33 <- norMmix(
    name = "#3 3d EEI",
    weight = c(0.3, 0.4, 0.3),
    mu = matrix(c(0,0,0,2,0,0,5,0,0),3,3),
    Sigma = array(rep(c(3,0,0,0,1,0,0,0,2),3), c(3,3,3)),
    model = "EEI"
    )

MW34 <- norMmix(
    name = "#4 3d VEI",
    weight = c(0.1, 0.9),
    mu = matrix(rep(0,6), 3,2),
    Sigma = array(c(diag(1:3), 0.2*diag(3:1)), c(3,3,2)),
    model = "VEI"
    )


####----------------------------------------------------------------------------
## dim 5
####----------------------------------------------------------------------------

MW51 <- norMmix(
    name = "#1 5d EII",
    weight = 1:3/sum(1:3),
    mu = matrix(1:15, 5,3),
    Sigma = c(1),
    model = "EII"
    )




n2m <- function(obj) nMm2par(obj , trafo="clr1", obj$model)
ln2m <- function(obj) nMm2par(obj , trafo="logit", obj$model)
