context("test results of llnorMmix function for calculable cases")


test_that("EII test", {

          tr <- "clr1"
          m <- "EII"

          par. <- c(0,0,0)
          x <- rbind(c(0,0), c(1,1))
          p <- 2
          k <- 1

          ret <- llnorMmix(par., t(x), k, model=m)
          expect_equal(ret, (log(1/(2*pi))+log(1/(2*pi*exp(1)))))

          par. <- c(0,0,log(2))
          x <- rbind(c(0,0), c(1,1))
          p <- 2
          k <- 1

          ret <- llnorMmix(par., x, k, model=m)
          expect_equal(ret, 2*(-log(2*pi)- 0.5*log(4)) -0.5*1)

})


test_that("VII test", {

          tr <- "clr1"
          m <- "VII"

          par. <- c(0,0,0,1,1,0,log(2))
          x <- rbind(c(0,0), c(1,1), c(-1,-1))
          p <- 2
          k <- 2

          ret <- llnorMmix(par., t(x), k, model=m)
          expval <- 0.5*exp( -0.5*p*(log(2*pi)) -0.5*c(0,2,2)) + 0.5*exp( -0.5*p*(log(2*pi)+ log(2)) -0.5*c(1,0,4) )
          expect_equal(ret, sum(log(expval)))

})
