context("test llnorMmix against llmvtnorm to check values of llnorMmix")



test_that("EII llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EII"

          par. <- nc2p(MW21)
          x <- rnorMmix(511, obj=MW21)
          p <- 2
          k <- 1

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})



test_that("VII llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VII"

          set.seed(2019)
          par. <- nc2p(MW23)
          x <- rnorMmix(511, obj=MW23)
          p <- 2
          k <- 3

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})




test_that("EEI llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EEI"

          par. <- nc2p(MW26)
          x <- rnorMmix(511, obj=MW26)
          p <- 2
          k <- 2

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})


test_that("VEI llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VEI"

          par. <- nc2p(MW27)
          x <- rnorMmix(511, obj=MW27)
          p <- 2
          k <- 2

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})


test_that("EVI llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EVI"

          par. <- nc2p(MW28)
          x <- rnorMmix(511, obj=MW28)
          p <- 2
          k <- 3

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})


test_that("VVI llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VVI"

          par. <- nc2p(MW29)
          x <- rnorMmix(511, obj=MW29)
          p <- 2
          k <- 2

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})


test_that("EEE llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EEE"

          par. <- nc2p(MW210)
          x <- rnorMmix(511, obj=MW210)
          p <- 2
          k <- 2

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})


test_that("VEE llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VEE"

          par. <- nc2p(MW211)
          x <- rnorMmix(511, obj=MW211)
          p <- 2
          k <- 2

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})


test_that("EVV llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "EVV"

          par. <- nc2p(MW212)
          x <- rnorMmix(511, obj=MW212)
          p <- 2
          k <- 2

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})


test_that("VVV llnorMmix against llmvtnorm", {
          tr <- "clr1"
          mo <- "VVV"

          par. <- nc2p(MW213)
          x <- rnorMmix(511, obj=MW213)
          p <- 2
          k <- 2

          retnMm <- llnorMmix(par., t(x), k, model=mo)
          retmvt <- llmvtnorm(par., x, k, model=mo)

          expect_equal(retnMm,retmvt)
})
