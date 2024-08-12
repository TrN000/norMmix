context("test norMmix class utilities")

expect_identical <- function(x, y) expect_true(identical(x, y))

test_that("we can convert models from `nor1mix` to `norMmix`.", {
  skip_if_not_installed("nor1mix")
  library(nor1mix)
  original <- nor1mix::MW.nm4
  coerced <- nor1toMmix(original)
  expect_equal(coerced$mu, matrix(0, 1, 2))
  expect_identical(coerced$dim, 1L)
  expect_true(inherits(coerced, "norMmix"))
})

test_that("constructor handles p=1 gracefully", {

    ## norMmix() might confuse cases where p=1 with k=1, mistaking a matrix for a vector.
    nMm <- norMmix(mu=matrix(c(1,1,2), 1, 3), # 3 means with dim 1
            model="EII")
    expect_identical(nMm$dim, 1L)
    expect_identical(nMm$k, 3L)

    nMm <- norMmix(mu=c(1,1,3), # 1 mean with dim 3, tacitly converted to matrix(mu, length(mu), 1)
                   model="EII")
    expect_identical(nMm$dim, 3L)
    expect_identical(nMm$k, 1L)
})


test_that("constructor builds and throws correctly for varying Sigmas.", {
    mu <- matrix(0,2,2) # means don't matter here
    sig <- 2

    ## should correctly build out of single value
    nMm <- norMmix(mu, Sigma=sig)
    expect_equal(max(abs(sig*diag(2) - nMm$Sigma[,,1])), 0)
    ## TODO: finish regression test
})
