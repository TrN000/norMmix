context("test input/output of par2nMm")



test_that("io test using MWnm", {
          
          m <- MW21
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)

          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW23
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)

          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW24
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW22
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW25
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW26
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW27
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW28
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW29
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW210
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW211
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW212
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)

          m <- MW213
          a <- nc2p(m)
          b <- par2nMm(a,m$dim,m$k,m$model)
          expect_equal(m$weight,b$weight)
          expect_equal(m$mu,b$mu)
          expect_equal(m$Sigma[,,1],b$Sigma[,,1])
          expect_equal(m$k,b$k)
          expect_equal(m$dim,b$dim)







})
