## functions that handle RDS files, that contain fittednorMmix objects

# extract BIC from files, return array with sensible dimnames
massbic <- function(string, DIR) {

    nm1 <- readRDS(file=file.path(DIR,string[1]))
    cl <- nm1$fit$k
    mo <- nm1$fit$models

    val <- array(0, lengths(list(cl, mo, string)))

    for (i in 1:length(string)) {
        nm <- readRDS(file=file.path(DIR,string[i]))
        val[,,i] <- BIC(nm$fit)[[1]]
    }
    dimnames(val) <- list(clusters=cl, models=mo, simulation=string)

    ## TODO: add more info to result, give it a class maybe
    val
}

# applies mclust along a model selection of fitnMm
massbicm <- function(string, DIR) {
    nm <- readRDS(file.path(DIR, string[1]))
    cl <- nm$fit$k
    mo <- nm$fit$models
    valm <- array(0, lengths(list(cl,mo,string)))
    for (i in 1:length(string)) {
        nm <- readRDS(file.path(DIR, string[i]))
        x <- nm$fit$x
        valm[,,i] <- mclust::Mclust(x, G=cl, modelNames=mo)$BIC
    }
    dimnames(valm) <- list(clusters=cl, models=mo, files=string)
    -valm
}

# plot array of massbic
massplot <- function(f, main="unnamed", p) {
    ran <- extendrange(f)
    size <- dim(f)[3]
    cl <- as.numeric(dimnames(f)$clusters)
    adj <- exp(-0.002*size) # set so at n=1000 alpha value is ~0.1
    models <- mods()
    ## FIXME: should not assume all models present.
    op <- sfsmisc::mult.fig(mfrow=c(4,5), main=main, mar= 0.1 +c(2,4,4,1))
    for (i in 1:10) {
        if (!missing(p)) {
            matplot(f[,i,], lty=1, col=adjustcolor(rainbow(10)[i],adj), 
                    type="l", ylim=ran,main=models[i])
            axis(3, at=seq_along(cl), 
                 label=npar(cl,p,models[i]))
        } else {
            matplot(f[,i,], lty=1, col=adjustcolor(rainbow(10)[i],adj), 
                    main=models[i], type="l", ylim=ran)
        }
    }
    for (i in 1:10) {
        boxplot(t(f[,i,]), lty=1, col=adjustcolor(rainbow(10)[i],0.4), 
                main=models[i], type="l", ylim=ran)
    }
    par(op$old.par)
}

# compare two massbic arrays
compplot <- function(f, g) {
    ylim <- extendrange(c(f,g))
    adj <- 0.4
    op <- par(mfrow=c(2,5))
    for (i in 1:10) {
        matplot(f[,i,], lty=1, col=adjustcolor(rainbow(20)[i],adj), 
                main=models[i], type="l", ylim=ylim)
        matplot(g[,i,], lty=1, col=adjustcolor(rainbow(20)[i+10],adj), 
                main=models[i], type="l", ylim=ylim, add=TRUE)
    }
    par(op)
}
