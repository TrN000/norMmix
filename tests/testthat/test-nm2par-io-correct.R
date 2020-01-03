context("test io correctness of nMm2par")


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("test length", {
          t1 <- nMm2par(MW21, model=MW21$model)
          expect_equal(length(t1),3)
          t2 <- nMm2par(MW23, model=MW23$model)
          expect_equal(length(t2),11)
          t4 <- nMm2par(MW24, model=MW24$model)
          expect_equal(length(t4),7)
          t6 <- nMm2par(MW26, model=MW26$model)
          expect_equal(length(t6),7)
          t7 <- nMm2par(MW27, model=MW27$model)
          expect_equal(length(t7),8)
          t8 <- nMm2par(MW28, model=MW28$model)
          expect_equal(length(t8),12)
          t9 <- nMm2par(MW29, model=MW29$model)
          expect_equal(length(t9),9)
          t10 <- nMm2par(MW210, model=MW210$model)
          expect_equal(length(t10),8)
          t11 <- nMm2par(MW211, model=MW211$model)
          expect_equal(length(t11),9)
          t12 <- nMm2par(MW212, model=MW212$model)
          expect_equal(length(t12),10)
          t13 <- nMm2par(MW213, model=MW213$model)
          expect_equal(length(t13),11)
          
})




test_that("test if nMm2par(p2n()) == id ", {

          tr <- "clr1"
          m <- MW211
          l <- nMm2par(MW211, MW211$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW212
          l <- nMm2par(MW212, MW212$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW213
          l <- nMm2par(MW213, MW213$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW26
          l <- nMm2par(MW26, MW26$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW27
          l <- nMm2par(MW27, MW27$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW28
          l <- nMm2par(MW28, MW28$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW29
          l <- nMm2par(MW29, MW29$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW21
          l <- nMm2par(MW21, MW21$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW22
          l <- nMm2par(MW22, MW22$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW23
          l <- nMm2par(MW23, MW23$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW24
          l <- nMm2par(MW24, MW24$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

          m <- MW25
          l <- nMm2par(MW25, MW25$model)
          k <- nMm2par(par2nMm(l,m$dim,m$k,m$model),
                   model=m$model)
          expect_equal(k,l)

})




test_that("test wrong inputs in mu", {
          
          tr <- "clr1"
          o1 <- MW210
          o1$mu <- cbind( c(0,0),c(0,0),c(0,0) )
          expect_error(nMm2par(o1, model=o1$model))

          o1$mu <- cbind( c(0,0,0),c(0,0,0) )
          expect_error(nMm2par(o1, model=o1$model))

          o1$mu <- cbind( c(0,0),c(0,"hello") )
          expect_error(nMm2par(o1, model=o1$model))

          o1$mu <- cbind( c(0,0),c(0,NA) )
          expect_error(nMm2par(o1, model=o1$model))
})


test_that("test wrong inputs in w", {
          tr <- "clr1"
          o2 <- MW28

          o2$weight <- c(1,1)
          expect_error(nMm2par(o2, model=o1$model))
          
          o2$weight <- c(1,1)/2
          expect_error(nMm2par(o2, model=o1$model))
          
          o2$weight <- c(1,1,1,1)/4
          expect_error(nMm2par(o2, model=o1$model))



})




