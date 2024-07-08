context("test that assertion function `sigmaAgainstModel` correctly throws on incorrect Sigma")

# Note: we do not assert on non-sym. pos. def. as that check already
# happens in okSigma.


# expect_message
#
# Our function has return type Either(Bool, [Char]) so we want to test
# that the right message is returned on error.  This helper function
# does just that.  In the boolean case, we can use expect_true.
expect_message <- function(ret, msg) {
    expect_type(ret, "character")
    if (ret != TRUE) {
        expect_true(startsWith(ret, msg), info=ret)
    }    
}

## EII
#
# EII means constant volume \alpha across all cov. mat. and the shape
# \Lambda is an identity matrix.  we test for equal volume, constant
# diagonal and no off-diagonal coeffs.

test_that("EII throws on off diagonal components.", {
    sig <- array(c(diag(c(2,2))),  c(2,2,2))
    sig[1,2,1] <- 2
    ret <- sigmaAgainstModel(sig, "EII")
    expect_message(ret, "Sigma has off-diagonal coefficients")
})

test_that("EII throws when diagonal entries are not equal.", {
    sig <- array(c(diag(c(1,1,2))),  c(3,3,4))
    ret <- sigmaAgainstModel(sig, "EII")
    expect_message(ret, "Sigma is not EII")
})

test_that("EII accepts correctly", {
    sig <- array(c(3*diag(3)), c(3,3,4))
    ret <- sigmaAgainstModel(sig, "EII")
    expect_true(ret)            
})


## VII tests
#
# VII is the same as EII, except the volumes may vary.  So we only
# need to verify, that the diagonals are constant and not that they
# are all equal.

test_that("VII throws when diagonal entries are not equal.", {
    sig <- array(c(diag(c(1,1,2))),  c(3,3,4))
    ret <- sigmaAgainstModel(sig, "VII")
    expect_message(ret, "Sigma is not VII")
})


test_that("VII throws on nonzero off diagonals.", {
    sig <- array(3*diag(2), c(2,2,2))
    sig[1,2,1] <- 42
    ret <- sigmaAgainstModel(sig, "VII")
    expect_message(ret, "Sigma has off-diagonal coefficients")
})

test_that("VII accepts correct model", {
    sig <- array(3*diag(2), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "VII")
    expect_true(ret)
})


# EEI
#
# EEI deviates from the previous models, in that now a varying
# diagonal entries are permitted.  The model separates volume and
# shape into the variables \alpha and \Lambda, so \Lambda is required
# to have det=1.  In this case we check if all diagonals are equal and
# that correctly asserts.

test_that("EII throws on nonzero off diagonals", {
    sig <- array( diag(c(1,2)), c(2,2,2))
    sig[1,2,1] <- 1
    ret <- sigmaAgainstModel(sig, "EEI")
    expect_message(ret, "Sigma has off-diagonal coefficients")
})

test_that("EII throws on unequal cov mat", {
    sig <- array( c(
        diag(c(1,2)),
        diag(c(2,3))), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EEI")
    expect_message(ret, "Sigma is not EEI")
})

test_that("EEI accepts correct model", {
    sig <- array( diag(c(1,2)), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EEI")
    expect_true(ret)
})


# VEI

test_that("VEI throws on off-diagonal entries", {
    sig <- array( c(
        2*diag(c(2,3)),
        2*diag(c(2,3))), c(2,2,2))
    sig[1,2,1] <- 1
    ret <- sigmaAgainstModel(sig, "VEI")
    expect_message(ret, "Sigma has off-diagonal coefficients")
})

test_that("VEI throws when covariances are not multiple of each other.", {
    sig <- array( c(
        2*diag(c(2,3)),
        32*diag(c(2,3))), c(2,2,2))
    sig[2,2,2] <- 97
    ret <- sigmaAgainstModel(sig, "VEI")
    expect_message(ret, "Sigma is not VEI. Matrix at")
})

test_that("VEI accepts correct model", {
    sig <- array( c(
        2*diag(c(2,3)),
        32*diag(c(2,3))), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "VEI")
    expect_true(ret)
})


## EVI

test_that("EVI throws on off-diagonal components", {
    sig <- array( c(
        diag(c(1,1)),
        diag(c(2,3))), c(2,2,2))
    sig[1,2,1] <- 3
    ret <- sigmaAgainstModel(sig, "EVI")
    expect_message(ret, "Sigma has off-diagonal coefficients")
})

test_that("EVI throws on non-equal volume", {
    sig <- array( c(
        diag(c(1,2)),
        diag(c(2,4))), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EVI")
    expect_message(ret, "Sigma does not have equal volume")
})

test_that("EVI accepts correct model", { 
    sig <- array( c(
        diag(c(1,1)),
        diag(c(2,0.5))), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EVI")
    expect_true(ret)
})


## VVI
#
# Relatively easy to test. This one just needs to have no off diagonal
# components.

test_that("VVI throws on non-zero off-diagonals", {
    sig <- array(c(
        diag(2,2),
        diag(3,2)), c(2,2,2))
    sig[1,2,1] <- 1
    ret <- sigmaAgainstModel(sig, "VVI")
    expect_message(ret, "Sigma has off-diagonal coefficients")
})

test_that("VVI accepts correct model", {
    sig <- array(c(
        diag(2,3),
        diag(3,4)), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "VVI")
    expect_true(ret)
})


#####################
# non-diagonal models
#####################


## EEE
#
# verify that all cov matrices are equal

test_that("EEE throws on unequal mat", {
    sig <- array(c(1,2,3,4, 1,2,3,5), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EEE")
    expect_message(ret, "Sigma not equal")
})

test_that("EEE accepts correct model", {
    sig <- array(c(1,2,3,4), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EEE")
    expect_true(ret)
})


## VEE
#
# verify that all cov matrices have the same shape

test_that("VEE throws on unequal shape", {
    sig <- array(c(
        c(1,2,3,4),
        2*c(1,2,3,5)
    ), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "VEE")
    expect_message(ret, "")
})

test_that("VEE accepts correct model", {
    sig <- array(c(
        c(1,2,3,4),
        2*c(1,2,3,4)
    ), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "VEE")
    expect_true(ret)
})


## EVV
#
# verify that all cov mat have the same volume

test_that("EVV throws on unequal volume", {
    sig <- array(c(
        c(1,2,3,4),
        2*c(1,3,2,4)
    ), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EVV")
    expect_message(ret, "")
})

test_that("EVV accepts correct model", {
    sig <- array(c(
        c(1,2,3,4),
        c(1,3,2,4)
    ), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EVV")
    expect_true(ret)
})

## VVV
#
# don't verify anything
test_that("VVV accepts any model", {
    sig <- array(c(1,2,3,4,  6,5,4,3), c(2,2,2))
    ret <- sigmaAgainstModel(sig, "VVV")
    expect_true(ret)
})
