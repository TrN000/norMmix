## functions that handle RDS files, that contain fittednorMmix objects


# extract BIC from files, return array with sensible dimnames
massbic <- function(string, DIR) {
    nm1 <- readRDS(file=file.path(DIR,string[1]))
    cl <- nm1$fit$k
    mo <- nm1$fit$models

    val <- array(0, lengths(list(cl, mo, string)))
    dims <- vector(mode="integer", length=length(string))

    for (i in 1:length(string)) {
        nm <- readRDS(file=file.path(DIR,string[i]))
        val[,,i] <- BIC(nm$fit)[[1]]
        dims[i] <- nm$fit$p
    }
    dimnames(val) <- list(components=cl, models=mo, simulation=string)
    attr(val, "dims") <- dims

    ## TODO: add more info to result, give it a class maybe
    val
}


# applies mclust along a model selection of fitnMm
massbicm <- function(string, DIR) {
    nm <- readRDS(file.path(DIR, string[1]))
    x1 <- nm$fit$x
    cl <- nm$fit$k
    mo <- nm$fit$models
    valm <- array(0, lengths(list(cl,mo,string)))
    dims <- vector(mode="integer", length=length(string))
    for (i in 1:length(string)) {
        nm <- readRDS(file.path(DIR, string[i]))
        x <- nm$fit$x
        stopifnot(all.equal(x,x1))
        valm[,,i] <- mclust::Mclust(x, G=cl, modelNames=mo)$BIC
        dims[i] <- nm$fit$p
    }
    dimnames(valm) <- list(components=cl, models=mo, files=string)
    attr(valm, "dims") <- dims
    -valm
}


# plot array of massbic
massplot <- function(f, main="unnamed", 
                     adj=exp(-0.002*size), col=nMmcols[1],
                     mar=0.1+c(1.4,2,3,1),
                  ...) {
    ran <- extendrange(f)
    size <- dim(f)[3]
    cl <- as.numeric(dimnames(f)$components)
    p <- attr(f, "dims")
    adj <- adj 
    models <- mods()
    ## FIXME: should not assume all models present.
    op <- sfsmisc::mult.fig(mfrow=c(4,5), main=main, mar=mar)
    for (i in 1:10) {
        if (!is.null(p)) {
            matplot(f[,i,], lty=1, col=adjustcolor(col,adj), 
                    type="l", ylim=ran, ylab='', ...)
            title(main=models[i], line=2)
            axis(3, at=seq_along(cl), 
                 labels=dfnMm(cl,p[1],models[i]))
        } else {
            matplot(f[,i,], lty=1, col=adjustcolor(col,adj), 
                    type="l", ylim=ran, ylab='',  ...)
            title(main=models[i], line=3)
        }
    }
    for (i in 1:10) {
        boxplot(t(f[,i,]), lty=1, col=adjustcolor(col,0.4), 
                main=models[i], type="l", ylim=ran)
    }
    par(op$old.par)
}


# compare two massbic arrays
compplot <- function(f, g, h=NULL, main="unnamed", 
                     adj=1/dim(f)[3], col=nMmcols[1:3], 
                     mar=0.1+c(1.4,2,3,1),
                     compnames = c("clara", "mclVVV", "Mclust"),
                     oma=c(7,0,2.8,0),
                  ...) {
    ylim <- extendrange(c(f,g))
    adj <- 0.4
    if (is.null(compnames)) {oma <- c(1,0,3,0)}
    op <- sfsmisc::mult.fig(mfrow=c(2,5), main=main, mar=mar, oma=oma)
    models <- dimnames(f)$models
    cl <- as.numeric(dimnames(f)$components)
    p <- attr(f, "dims")
    for (i in 1:10) {
        matplot(f[,i,], lty=1, col=adjustcolor(col[1],adj), 
                main=models[i], type="l", ylim=ylim, ylab='', ...)
        axis(3, at=seq_along(cl), 
             labels=dfnMm(cl,p[1],models[i]))
        matplot(g[,i,], lty=1, col=adjustcolor(col[2],adj), 
                main=models[i], type="l", ylim=ylim, add=TRUE, ylab='', ...)
        if (!is.null(h)) {
            matplot(h[,i,], lty=1, col=adjustcolor(col[3],adj), 
                    main=models[i], type="l", ylim=ylim, add=TRUE, ylab='', ...)
        }
    }
    par(op$old.par)
    if (!is.null(compnames)) {
        op <- par(xpd=TRUE)
        legend("bottom", horiz=TRUE, legend=compnames, fill=col, inset=c(0, -0.05))
        par(op)
    }
}


epfl <- function(files, savdir, subt=11, ...) {
    # files: list of character vectors
    # savdir: string, specifies directory
    # subt: how many characters to subtract from filename
    #       for 'main'
    
    stopifnot(is.list(files), dir.exists(savdir))
    setwd(savdir)
    for (fi in files) {
        ## plot
        if (length(fi)==0) {next}
        main <- substring(fi[1], 1, nchar(fi[1])-subt)
        f <- massbic(fi, savdir)
        g <- massbicm(fi, savdir)
        saveRDS(g, file=file.path(savdir, paste0(main, "_mcl.rds")))
        pdf(file=paste0(main,".pdf"))
        massplot(f, main=main, ...)
        dev.off()
        pdf(file=paste0(main,"_mcl.pdf"))
        massplot(g, main=paste0(main, "_mcl"), ...)
        dev.off()
        #pdf(file=paste0(main,"_comp.pdf"))
        #compplot(f,g, main=paste0(main, "_comp"), ...)
        #dev.off()
    }
}
