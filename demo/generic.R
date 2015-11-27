old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

if( !require('ggplot2') ) {
  stop('This demo requires package ggplot2. Please install.')
}

#### Generic component ####

## The generic component admits custom incidence and covariance or precision
## matrices. It can be used to manually incorporate features not included.
## We can check that it gives the same results than existing components.


data(globulus)

## LMM
res.lmm <- remlf90(fixed  = phe_X ~ gg + x,
                   random = ~ bl,
                   data   = globulus)

bl.inc <- model.matrix(~ 0 + bl, globulus)
bl.cov <- diag(nlevels(globulus$bl))
res.lmm2 <- remlf90(fixed   = phe_X ~ gg + x,
                    generic = list(bl = list(bl.inc,
                                             bl.cov)),
                    data    = globulus)

all.equal(res.lmm, res.lmm2, check.attributes = FALSE)


## Animal model
ped <- build_pedigree(c('self', 'dad', 'mum'),
                      data = globulus)

system.time(
  res.anm <- remlf90(fixed  = phe_X ~ gg + x,
                     random = ~ bl,
                     genetic = list(model    = 'add_animal',
                                    pedigree = ped,
                                    id       = 'self'),
                     data   = globulus)
)  ## ~ 0.16


gen.inc <- model.matrix(res.anm)$genetic
gen.cov <- pedigreemm::getA(ped)
system.time(
  res.anm2 <- remlf90(fixed   = phe_X ~ gg + x,
                      random = ~ bl,
                      generic = list(gen = list(gen.inc,
                                                gen.cov)),
                      data    = globulus)
)  ## ~ 0.48

all.equal(res.anm, res.anm2, check.attributes = FALSE, tol = 1e-07)


## generic admits a list of random effects
system.time(
  res.anm3 <- remlf90(fixed   = phe_X ~ gg + x,
                    generic = list(bl = list(bl.inc,
                                             bl.cov),
                                   gen = list(gen.inc,
                                              gen.cov)),
                    data    = globulus)
)  ## ~ 0.625

all.equal(res.anm, res.anm3, check.attributes = FALSE, tol = 1e-07)

## All give the same results, but more specific
## effects are more efficiently implemented


## Spatial models

# The genetic component
gen.globulus <- list(model    = 'add_animal', 
                     pedigree = globulus[,1:3],
                     id       = 'self')

system.time(
  res.spl  <- remlf90(fixed  = phe_X ~ gg,
                      genetic = gen.globulus, 
                      spatial = list(model = 'splines', 
                                     coord = globulus[, c('x','y')], 
                                     n.knots = c(7, 7)), 
                      data = globulus,
                      method = 'em')
)  # ~ 6 s

spl.inc <- model.matrix(res.spl)$spatial
spl.cov <- get_structure(res.spl)$spatial

system.time(
  res.spl2  <- remlf90(fixed  = phe_X ~ gg,
                       genetic = gen.globulus, 
                       generic = list(spl = list(spl.inc,
                                                 spl.cov)), 
                       data = globulus,
                       method = 'em')
)  # ~ 6 s

all.equal(res.spl, res.spl2, check.attributes = FALSE)

system.time(
  res.ar  <- remlf90(fixed  = phe_X ~ gg,
                     genetic = gen.globulus, 
                     spatial = list(model = 'AR', 
                                    coord = globulus[, c('x','y')],
                                    rho = c(.85, .8)), 
                     data = globulus)
)  # ~ .6


ar.inc <- model.matrix(res.ar)$spatial
ar.prc <- get_structure(res.ar)$spatial

## Note that the latter is a precision matrix (inverse covariance)
## So we need to tell so explicitly:
system.time(
  res.ar2  <- remlf90(fixed  = phe_X ~ gg,
                      genetic = gen.globulus, 
                      generic = list(ar = list(ar.inc,
                                               prec = ar.prc)), 
                      data = globulus)
)  # ~ .96

all.equal(res.ar, res.ar2, check.attributes = FALSE)


## Being able to include several independent generic effects
## lets you fit complementary spatial models
## (for example, to capture long and short-range spatial effects)

res.sp2 <- remlf90(fixed  = phe_X ~ gg,
                   genetic = gen.globulus, 
                   spatial = list(model = 'AR', 
                                  coord = globulus[, c('x','y')],
                                  rho = c(.3, .3)), 
                   generic = list(spl = list(spl.inc,
                                             spl.cov)),
                   method  = 'em',
                   data = globulus)


## compare the variograms of residuals
variogram(res.spl, 'isotropic')
variogram(res.sp2, 'isotropic')

## the AR component in the combined model captures the short-range
## autocorrelation still present in the residuals from the splines-only model
compare.plots(list(
  'Residuals splines-only model' = plot(res.spl, 'resid'),
  'short-range AR from combined model' = plot(res.sp2, 'spatial')))

