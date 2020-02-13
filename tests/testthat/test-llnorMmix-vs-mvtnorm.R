context("test llnorMmix against llmvtnorm to check values of llnorMmix")


test_that("EII llnorMmix against llmvtnorm", {
          mo <- "EII"
          par. <- nMm2par(MW21, model=mo)
          x <- rnorMmix(511, obj=MW21)
          k <- MW21$ k # 1
          ##
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          ##
          expect_equal(retnMm,retmvt)
})


test_that("VII llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VII"
          set.seed(2019)
          par. <- nMm2par(MW23, model=mo)
          x <- rnorMmix(511, obj=MW23)
          k <- MW23$ k
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})


test_that("EEI llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EEI"
          par. <- nMm2par(MW26, model=mo)
          x <- rnorMmix(511, obj=MW26)
          k <- MW26$ k
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})


test_that("VEI llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VEI"
          par. <- nMm2par(MW27, model=mo)
          x <- rnorMmix(511, obj=MW27)
          k <- 2
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})


test_that("EVI llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EVI"
          par. <- nMm2par(MW28, model=mo)
          x <- rnorMmix(511, obj=MW28)
          k <- 3
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})


test_that("VVI llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VVI"
          par. <- nMm2par(MW29, model=mo)
          x <- rnorMmix(511, obj=MW29)
          k <- 2
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})


test_that("EEE llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EEE"
          par. <- nMm2par(MW210, model=mo)
          x <- rnorMmix(511, obj=MW210)
          k <- 2
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})


test_that("VEE llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VEE"
          par. <- nMm2par(MW211, model=mo)
          x <- rnorMmix(511, obj=MW211)
          k <- 2
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})


test_that("EVV llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EVV"
          par. <- nMm2par(MW212, model=mo)
          x <- rnorMmix(511, obj=MW212)
          k <- 2
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})


test_that("VVV llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VVV"
          par. <- nMm2par(MW213, model=mo)
          x <- rnorMmix(511, obj=MW213)
          k <- 2
          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par.,   x , k, model=mo)
          expect_equal(retnMm,retmvt)
})
