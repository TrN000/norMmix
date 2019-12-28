norMmix package
=====

This package implements a model fitting algorithm for normal mixture models.

TODO:
- [X] trim the fat; remove simulation functionality
- [X] make init method optional with no default
- [] spurious cluster suppression
- [] add vignette
- [] test against windows
- [] maybe add gitwiki tutorial
- [] change norMmixMLE to class c("norMmixMLE", "norMmix")
- []


file structure:
        R
        ├── Cholesky.R          implements ldlt decomp
        ├── llnorMmix.R         log-likelihood function
        ├── norMmixMLE.R        MLE algorithm
        ├── norMmix.R           base S3 class + some convenience functions
        ├── param.R             parametrization functions
        ├── plot.R              plot methods
        └── zmarrwandnMm.R      example objects
