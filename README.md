norMmix package
=====

This package implements a model fitting algorithm for normal mixture models.

TODO:
- [X] make init method optional with no default
- [ ] spurious cluster suppression
- [X] (WIP) add vignette
- [X] test against windows
- [X] (WIP) reintroduce functionality for many MLE passes for model selection


## Install

On linux machines, using R only:

```{bash}
    git clone https://github.com/TrN000/norMmix.git
    R CMD build norMmix
    R CMD INSTALL norMmix_0.0-2.tar.gz
```

To build and check the package the following packages are necessary.

```{R}
    > install.packages(
        c("mvtnorm", "mclust", "sfsmisc", "nor1mix", "testthat",
          "knitr", "rmarkdown", "compositions", "mixtools"))
```
