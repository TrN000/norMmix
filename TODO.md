

0. `ssClara2kL()` {was 'ssClaraL'}: seems to give too small samples in
	several of the BSc thesis simulations.  Provide "better" alternative!

1. Provide  as.norMmix() generic and   as.norMmix.nor1mix()  to get 1-D
   comparisons and MW<n>  from CRAN package  `nor1mix`

3. `norMmix()` has trouble with 1 dimensional and/or 1 component mixtures.
   deciding criteria needed.

4. rename nMmcol to Trubetskoy10

6. Documentation for plot methods now completely broken. Arguments like
   newWindow no longer supported. Needs to be rewritten for 2d and >2d.

7. norMmix(...) creation does *NOT* check `Sigma` in case of restricted
   parametrizations,  but just assumes the caller of norMmix() makes no
   mistake.  This is clearly too optimistic!
   -> split into "private" and public method, where the public method
   does argument checking and private one assumes correct args.
   Should use the private one in MLE algorithm as it gets called there in a
   loop.

8. Currently always use *full* covariance parameter Sigma, even in cases
   such as EII etc.  Should we allow *both* (full and minimal) parametrizations

9. `manyMLE()` is not yet exported & documented; before doing so, needs
   tweaks (by MM):
	- `savdir` and `name` with default name to save as `"*rds"`: nice idea but
      should *not* be part of `manyMLE` but separate small utility.
	- modify argument checking; use namespace global
	  `norModels <- eval(formals(norMmix)$model)`

10. *BUG* in either `npar()` or `nMm2par()` -- see 'FIXME' in `man/nMm2par.Rd`

==============================================================================
DONE:

2. Also port most of the "fit.R" from `norMmix_Bthesis`; i.e., the  `fitnMm()`
  function {with better name!} -->
  `~/Betreute-Arbeiten/NicolasTrutmann/BSc_thesis+MM/norMixBthesis/R/fit.R`
  and then all  the   `<fun>.fittednorMmix()`  methods
  {but probably also change the name of the class from
  `"fittednorMmix"` to something like
  `"manyNormixMLE"`
     --> Nicolas Trutmann has done it (Jun 21 2020): name  `manyMLE()` ==> `~/R/D/GH/norMmix/R/fit.R`

4. norMmix.Rd does not document use of non-array covar. mats. as init. values

5. decide on par args in ndplot and how to leave it exposed to the user.
   maybe do.call(par, parargs) construct?? and put parargs=NULL in arguments.
   user can then overwrite first call to par().

   ==> MM: Using sfsmisc::mult.fig()  and allow the user to change defaults
   for  mult.fig() is easier.

  ===> **MUCH better** is really to change the setup completely, use
	   graphics :: pairs.default(.)
    by providing a correct  panel = function(.)    !!!

0. --> `R/plot.R`: FIXME: plot2d() <--> plotnd() are *NOT* compatible in their defaults

