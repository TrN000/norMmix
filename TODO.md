

1. Provide  as.norMmix() generic and   as.norMmix.nor1mix()  to get 1-D
   comparisons and MW<n>  from CRAN package  `nor1mix`

2. Also port most of the "fit.R" from `norMmix_Bthesis`; i.e., the  `fitnMm()`
  function {with better name!} -->
  `~/Betreute-Arbeiten/NicolasTrutmann/BSc_thesis+MM/norMixBthesis/R/fit.R`
  and then all  the   `<fun>.fittednorMmix()`  methods
  {but probably also change the name of the class from
  `"fittednorMmix"` to something like
  `"manyNormixMLE"`
