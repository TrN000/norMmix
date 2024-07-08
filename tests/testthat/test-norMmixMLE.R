## regression tests for norMmixMLE.

context("regression tests for norMmixMLE")


test_that("a small example works as intended", {
    k <- 2
    sampleSize <- 100
    model <- "EII"
    nMm  <- norMmix(
       matrix(c(-3,0,3,0), 2, 2), ## simple 2 mixture model
       Sigma = NULL,
       weight = rep(1/k, k),
       name = NULL,
       model = model
    )

    x <- rnorMmix(100, nMm)

    mleResult <- norMmixMLE(x, model=model, k=k, maxit=1, optREPORT=1e6, trace=0)
    expect_equal(mleResult$npar, 6) # proper number of params
    expect_equal(mleResult$optr$convergence, 1) # should not have converged after 1 iter
    expect_equal(inherits(mleResult, "norMmixMLE"), TRUE) # return value has the right class
})


test_that("k=1 also works without confusing norMmix.", {

    sampleSize  <- 10
    model  <- "EII"

    nMm  <- norMmix(
        matrix(c(0,0), 2, 1),
        model = model)

    ## assert on number of components and dimensions
    expect_equal(nMm$k, 1)
    expect_equal(nMm$dim, 2)
})

