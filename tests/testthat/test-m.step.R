context("test m.step correctness")


test_that("test against mclust m.step", {

    set.seed(2019); x <- matrix(runif(20), 10, 2)

    clu <- cluster::clara(x,2)$clustering

    library(mclust)

    index <- matrix(0, 10,2)
    index[cbind(1:10,clu)] <- 1

    mc <- mclust::mstep("VVV", x, index)

    nm <- mstep.nMm(x, index)


    expect_equal(mc$parameters$pro, nm$weight)
    expect_equal(c(mc$parameters$mean), c(nm$mu))
    expect_equal(c(mc$parameters$variance$sigma), c(nm$Sigma))

})


test_that("test over MW* mcl against nMm", {

    MWdat <- Filter(function(.) is.norMmix(get(., "package:norMmix")),
                     ls("package:norMmix", pattern = "^MW[1-9]"))

    for (i in MWdat) {
        nMm <- get(i, "package:norMmix")
        p <- nMm$dim

        for(j in 2:3) { # done for 1:5, but too long
            set.seed(2014+j); x <- rnorMmix(200, nMm)

            for (k in 1:2) { # done for 1:7
                clu <- cluster::clara(x,k)$clustering
                index <- matrix(0, 200,k)
                index[cbind(1:200,clu)] <- 1

                mc <- mclust::mstep("VVV", x, index)
                nm <- mstep.nMm(x, index)

                expect_equal(mc$parameters$pro, nm$weight)
                expect_equal(c(mc$parameters$mean), c(nm$mu))
                expect_equal(c(mc$parameters$variance$sigma), c(nm$Sigma))
            }
        }
    }
})

