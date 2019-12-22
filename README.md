norMmix package
=====

This package implements a model fitting algorithm for normal mixture models.

TODO:
- [] trim the fat; remove simulation functionality
- [] make init method optional with no default
- [] maybe add gitwiki tutorial
- [] spurious cluster suppression
- [] add vignette
- [] test against windows
- []
- []
- []


file structure:
    R
    ├── Cholesky.R          implements ldlt decomp
    ├── llnorMmix.R         log-likelihood function
    ├── norMmixMLE.R        MLE algorithm
    ├── norMmix.R           base S3 class + some convenience functions
    ├── param.R             parametrization functions
    ├── plot.R              plot methods
    ├── weight.R            centered log ratio implementation
    └── zmarrwandnMm.R      example objects
