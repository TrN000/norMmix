context("test that assertion function `sigmaAgainstModel` correctly throws on incorrect Sigma")

test_that("EII throws correctly", {
    sig <- array(0, c(2,2,2))
    ret <- sigmaAgainstModel(sig, "EII")

    expect_type(ret, "character") # equivalent to error message
})
