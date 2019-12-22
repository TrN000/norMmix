context("test weight transformation corrctness")


test_that("test clr1 against clr1inv", {
          w <- 1
          expect_equal(w,clr1inv(clr1(w))) # test behaviour when clr1(w) = NULL
          w1 <- c(1,2)/3
          expect_equal(w1,clr1inv(clr1(w1)))
          w2 <- 1:10000/sum(1:10000)
          expect_equal(w2,clr1inv(clr1(w2)))

})
