context("norMmix class utilities")

library(nor1mix)
test_that("we can convert models from `nor1mix` to `norMmix`.", {
  original <- nor1mix::MW.nm4
  coerced <- nor1toMmix(original)
  expect_equal(coerced$mu, matrix(0, 1, 2))
  expect_equal(coerced$dim, 1)
  expect_equal(class(coerced), "norMmix")
})
