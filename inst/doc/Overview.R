## ----setup, include = FALSE----------------------------------------------
library(breedR)
library(ggplot2)

## ----dataset, results = 'asis', echo = FALSE-----------------------------
knitr::kable(head(globulus))

## ---- echo = FALSE-------------------------------------------------------
cat(attr(globulus, 'comment'), sep ='\n')
str(globulus, give.attr = FALSE)

## ----provenance-test-----------------------------------------------------
res <- remlf90(fixed  = phe_X ~ 1,
               random = ~ gg,
               data   = globulus)

## ----mixed-model-varini, eval = FALSE------------------------------------
#  res <- remlf90(fixed = phe_X ~ 1,
#                 random = ~ gg,
#                 var.ini = list(gg = 2, resid = 10),
#                 data = globulus)

## ----summary-------------------------------------------------------------
summary(res)

## ----extraction-1--------------------------------------------------------
fixef(res)
ranef(res)

## ----extraction-2, fig.width = 3, fig.height = 3-------------------------
qplot(
      fitted(res),
      globulus$phe_X) +
  geom_abline(intercept = 0,
              slope = 1,
              col = 'darkgrey')
str(resid(res))
extractAIC(res)
logLik(res)


## ----random-spec---------------------------------------------------------
random = ~ gg + factor(mum)  # note that mum is numeric

## ----interaction-standard------------------------------------------------
random = ~ gg * factor(mum)

## ----interaction-workaround----------------------------------------------
dat <- transform(globulus,
                 interaction = factor(gg:bl))
random = ~ gg + bl + interaction

## ----hierarchical-exercise, purl=TRUE, message=FALSE---------------------
res.h <- remlf90(fixed = phe_X ~ 1,
                 random = ~ factor(mum) + gg,
                 data = globulus)


## ----factorial-exercise, purl=TRUE, message=FALSE------------------------
# Interaction variable
globulus.f <- transform(globulus,
                        gg_bl = factor(gg:bl))

res.f <- remlf90(fixed = phe_X ~ 1,
                 random = ~ gg + bl + gg_bl,
                 data = globulus.f)


## ----interaction-exercise, purl=TRUE-------------------------------------
summary(res)
summary(res.h)

## ----factorial-summary, purl=TRUE----------------------------------------
summary(res.f)

## ----AIC-factorial, message=FALSE, purl=TRUE-----------------------------
## result without interaction
res.f0 <- remlf90(fixed  = phe_X ~ 1,
                  random = ~ gg + bl,
                  data = globulus)
paste('AIC:', round(extractAIC(res.f0)),
      'logLik:', round(logLik(res.f0)))

## ----genetic, message = FALSE--------------------------------------------
res.animal <- remlf90(fixed  = phe_X ~ 1,
                      random = ~ gg,
                      genetic = list(model = 'add_animal', 
                                     pedigree = globulus[, 1:3],
                                     id = 'self'), 
                      data = globulus)

## ----genetic_result------------------------------------------------------
summary(res.animal)

## ----PBV, fig.show = 'hold'----------------------------------------------
## Predicted Breeding Values
# for the full pedigree first, and for the observed individuals
# by matrix multiplication with the incidence matrix
PBV.full <- ranef(res.animal)$genetic
PBV <- model.matrix(res.animal)$genetic %*% PBV.full

# Predicted genetic values vs.
# phenotype.
# Note: fitted = mu + PBV
qplot(fitted(res.animal), phe_X,
      data = globulus) +
  geom_abline(intercept = 0,
              slope = 1,
              col = 'gray')


## ----animal-residuals----------------------------------------------------
## Since coordinates have not
## been passed before they 
## must be provided explicitly.
coordinates(res.animal) <-
  globulus[, c('x', 'y')]
plot(res.animal, 'resid')

## ----genetic_variogram---------------------------------------------------
variogram(res.animal)

## ----variogram-isotropic, echo = FALSE-----------------------------------
variogram(res.animal, 
          coord = globulus[, c('x', 'y')],
          plot = 'isotropic')

## ----variogram-circle, echo = FALSE--------------------------------------
ang <- seq(0, 2*pi, by = pi/6)[- (1+3*(0:4))]
# Colours from
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")

coord <- data.frame(x = cos(ang),
                    y = sin(ang))
ggplot(coord, aes(x, y)) + 
  geom_point(size = 6, col = cbPalette[1]) + 
  geom_point(aes(x = 0, y = 0), size = 6) +
  coord_fixed() + 
  scale_x_continuous(breaks=NULL) + 
  scale_y_continuous(breaks=NULL) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


## ----variogram-heat, echo = FALSE----------------------------------------
variogram(res.animal, 
          coord = globulus[, c('x', 'y')],
          plot = 'heat')

## ----variogram-rect, echo = FALSE----------------------------------------
ggplot(coord, aes(x, y)) + 
  geom_point(size = 6, col = cbPalette[rep(c(1, 2, 2, 1), 2)]) + 
  geom_point(aes(x = 0, y = 0), size = 6) +
  coord_fixed() + 
  scale_x_continuous(breaks=NULL) + 
  scale_y_continuous(breaks=NULL) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


## ----variogram-anisotropic, echo = FALSE---------------------------------
variogram(res.animal, 
          coord = globulus[, c('x', 'y')],
          plot = 'aniso')

## ----variogram-direction, echo = FALSE-----------------------------------
ggplot(coord, aes(x, y)) + 
  geom_point(size = 6, col = cbPalette[rep(1:4, 2)]) + 
  geom_point(aes(x = 0, y = 0), size = 6) +
  coord_fixed() + 
  scale_x_continuous(breaks=NULL) + 
  scale_y_continuous(breaks=NULL) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

## ----spatial-blocks, message = FALSE-------------------------------------
# The genetic component (DRY)
gen.globulus <- list(model    = 'add_animal', 
                     pedigree = globulus[, 1:3],
                     id       = 'self')

res.blk <- remlf90(fixed   = phe_X ~ 1,
                   random  = ~ gg,
                   genetic = gen.globulus, 
                   spatial = list(model = 'blocks', 
                                  coord = globulus[, c('x', 'y')],
                                  id = 'bl'),
                   data    = globulus)

## ----genetic-spatial_result----------------------------------------------
summary(res.blk)

## ----genetic-spatial_variogram, fig.width = 6, fig.height = 5, echo = FALSE----
variogram(res.blk)

## ----splines, message = FALSE--------------------------------------------
## Use the `em` method! `ai` does not like splines
res.spl  <- remlf90(fixed   = phe_X ~ 1,
                    random  = ~ gg,
                    genetic = gen.globulus, 
                    spatial = list(model   = 'splines', 
                                   coord   = globulus[, c('x','y')]), 
                    data    = globulus, method  = 'em')

## ----globulus-fit, message = FALSE---------------------------------------
res.ar1  <- remlf90(fixed   = phe_X ~ 1,
                    random  = ~ gg,
                    genetic = gen.globulus, 
                    spatial = list(model = 'AR', 
                                   coord = globulus[, c('x','y')]), 
                    data    = globulus)


## ----change-residuals, fig.width = 8-------------------------------------
compare.plots(
  list(`Animal model only` = plot(res.animal, 'residuals'),
       `Animal/blocks model` = plot(res.blk, 'residuals'),
       `Animal/splines model` = plot(res.spl, 'residuals'),
       `Animal/AR1 model` = plot(res.ar1, 'residuals')))

## ----spatial-components, fig.width = 8-----------------------------------
compare.plots(list(Blocks  = plot(res.blk, type = 'spatial'),
                   Splines = plot(res.spl, type = 'spatial'),
                   AR1xAR1 = plot(res.ar1, type = 'spatial')))


## ----spatial-fullcomponents, fig.width = 8-------------------------------
compare.plots(list(Blocks  = plot(res.blk, type = 'fullspatial'),
                   Splines = plot(res.spl, type = 'fullspatial'),
                   AR1xAR1 = plot(res.ar1, type = 'fullspatial')))


## ----determine-nknots, fig.width = 5, fig.height = 2, echo = FALSE-------
ggplot(transform(data.frame(x = 10:1e3),
                 nok = breedR:::determine.n.knots(x)),
       aes(x, nok)) + 
  geom_line()

## ----spatial-exercise, eval = FALSE--------------------------------------
#  rho.grid <- expand.grid(rho_r = seq(.7, .95, length = 4),
#                          rho_c = seq(.7, .95, length = 4))

## ----spatial-exercise-1, purl = TRUE, message=FALSE----------------------
res.spl99  <- remlf90(fixed  = phe_X ~ 1, random = ~ gg,
                      genetic = gen.globulus,
                      spatial = list(model   = 'splines', 
                                     coord   = globulus[, c('x','y')],
                                     n.knots = c(9, 9)), 
                      data = globulus, method = 'em')

## ----spatial-exercise-1-results, purl = TRUE, echo = TRUE----------------
summary(res.spl)
summary(res.spl99)

## ----spatial-exercise-2, purl = TRUE-------------------------------------
qplot(rho_r, rho_c,
      fill = loglik,
      geom = 'tile',
      data = res.ar1$rho)

## ----spatial-exercise-3, purl = TRUE, message=FALSE----------------------
rho.grid <- expand.grid(rho_r = seq(.7, .95, length = 4),
                        rho_c = seq(.7, .95, length = 4))
res.ar.grid  <- remlf90(fixed  = phe_X ~ gg,
                        genetic = list(model = 'add_animal', 
                                       pedigree = globulus[,1:3],
                                       id = 'self'), 
                        spatial = list(model = 'AR', 
                                       coord = globulus[, c('x','y')],
                                       rho = rho.grid), 
                        data = globulus)
summary(res.ar.grid)

## ----competition-data----------------------------------------------------
# Simulation parameters
grid.size <- c(x=20, y=25) # cols/rows
coord <- expand.grid(sapply(grid.size,
                            seq))
Nobs <- prod(grid.size)
Nparents <- c(mum = 20, dad = 20)
sigma2_a <- 2   # direct add-gen var
sigma2_c <- 1   # compet add-gen var
rho      <- -.7 # gen corr dire-comp
sigma2_s <- 1   # spatial variance
sigma2_p <- .5  # pec variance
sigma2   <- .5  # residual variance

S <- matrix(c(sigma2_a,
              rho*sqrt(sigma2_a*sigma2_c),
              rho*sqrt(sigma2_a*sigma2_c),
              sigma2_c),
            2, 2)

set.seed(12345)
simdat <- 
  breedR.sample.phenotype(
    fixed   = c(beta = 10),
    genetic = list(model = 'competition',
                   Nparents = Nparents,
                   sigma2_a = S,
                   check.factorial=FALSE,
                   pec = sigma2_p),
    spatial = list(model = 'AR',
                   grid.size = grid.size,
                   rho   = c(.3, .8),
                   sigma2_s = sigma2_s),
    residual.variance = sigma2
    )

## Remove founders
dat <- subset(simdat,
              !(is.na(simdat$sire)
                & is.na(simdat$dam)))


## ----competition-fit, message = FALSE------------------------------------
system.time(
  res.comp <- remlf90(fixed   = phenotype ~ 1,
                      genetic = list(model = 'competition',
                                     pedigree = dat[, 1:3],
                                     id = 'self',
                                     coord = dat[, c('x', 'y')],
                                     competition_decay = 1,
                                     pec = list(present = TRUE)),
                      spatial = list(model = 'AR', 
                                     coord = dat[, c('x', 'y')],
                                     rho   = c(.3, .8)),
                      data = dat,
                      method = 'em')  # AI diverges
  )

## ----competition-results, results = 'asis', echo = FALSE, fig.width = 3----
var.comp <- 
  data.frame(True = c(sigma2_a, sigma2_c, rho, sigma2_s, sigma2_p, sigma2),
             Estimated = with(res.comp$var,
                              round(c(genetic[1, 1], genetic[2, 2],
                                      genetic[1, 2]/sqrt(prod(diag(genetic))),
                                      spatial, pec, Residual), digits = 2)),
             row.names = c('direct', 'compet.', 'correl.',
                          'spatial', 'pec', 'residual'))

knitr::kable(var.comp)

## ----competition-results-plot, echo = FALSE------------------------------
labels <- c(paste0(rep('sigma', 5),
                 c('[A]', '[C]', '[S]', '[P]', ''), '^2'),
            'rho')[c(1:2, 6, 3:5)]

ggplot(var.comp, aes(True, Estimated)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = 'darkgray') + 
  geom_text(label = labels, parse = TRUE, hjust = -.5) +
  expand_limits(x = 2.4, y = 2.6)


## ----competition-exercise-1, purl = TRUE, message=FALSE------------------

## compute the predicted effects for the observations
## by matrix multiplication of the incidence matrix and the BLUPs
pred <- list()
Zd <- model.matrix(res.comp)$'genetic_direct'
pred$direct <- Zd %*% ranef(res.comp)$'genetic_direct'

## Watch out! for the competition effects you need to use the incidence
## matrix of the direct genetic effect, to get their own value.
## Otherwise, you get the predicted effect of the neighbours on each
## individual.
pred$comp <- Zd %*% ranef(res.comp)$'genetic_competition'
pred$pec  <- model.matrix(res.comp)$pec %*% ranef(res.comp)$pec

## ----competition-exercise-1bis, purl = TRUE, echo = TRUE-----------------
comp.pred <-
  rbind(
    data.frame(
      Component = 'direct BV',
      True = dat$BV1,
      Predicted = pred$direct),
    data.frame(
      Component = 'competition BV',
      True = dat$BV2,
      Predicted = pred$comp),
    data.frame(
      Component = 'pec',
      True      = dat$pec,
      Predicted = as.vector(pred$pec)))

ggplot(comp.pred,
       aes(True, Predicted)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1,
              col = 'darkgray') +
  facet_grid(~ Component)

## ----generic-example, message = FALSE------------------------------------
## Fit a blocks effect using generic
inc.mat <- model.matrix(~ 0 + bl, globulus)
cov.mat <- diag(nlevels(globulus$bl))
res.blg <- remlf90(fixed  = phe_X ~ gg,
                   generic = list(block = list(inc.mat,
                                               cov.mat)),
                   data   = globulus)

## ----summary-generic, echo = FALSE---------------------------------------
summary(res.blg)

## ----prediction-remove-measure-------------------------------------------
rm.idx <- 8
rm.exp <- with(dat[rm.idx, ],
               phenotype - resid)
dat.loo <- dat
dat.loo[rm.idx, 'phenotype'] <- NA

## ----prediction-fit, echo = FALSE, message=FALSE-------------------------
res.comp.loo <- remlf90(fixed   = phenotype ~ 1,
                        genetic = list(model = 'competition',
                                       pedigree = dat[, 1:3],
                                       id = 'self',
                                       coord = dat[, c('x', 'y')],
                                       competition_decay = 1,
                                       pec = list(present = TRUE)),
                        spatial = list(model = 'AR', 
                                       coord = dat[, c('x', 'y')],
                                       rho   = c(.3, .8)),
                        data = dat.loo,
                        method = 'em')  

## ----prediction-validation, echo = FALSE, results='asis'-----------------
## compute the predicted effects for the observations
## by matrix multiplication of the incidence matrix and the BLUPs
Zd <- model.matrix(res.comp)$'genetic_direct'
pred.BV.loo.mat <- with(ranef(res.comp.loo), 
                        cbind(`genetic_direct`, `genetic_competition`))
pred.genetic.loo <- Zd[rm.idx, ] %*% pred.BV.loo.mat

valid.pred <- 
  data.frame(True = with(dat[rm.idx, ],
                         c(BV1, BV2, rm.exp)),
             Pred.loo = c(pred.genetic.loo,
                          fitted(res.comp.loo)[rm.idx]),
             row.names = c('direct BV', 'competition BV', 'exp. phenotype'))

knitr::kable(round(valid.pred, 2))

## ----prediction-exercise-1, purl = TRUE, results = 'asis', echo = -5-----
pred.BV.mat <- with(ranef(res.comp), 
                    cbind(`genetic_direct`, `genetic_competition`))

valid.pred$Pred.full <- c(Zd[rm.idx, ] %*% pred.BV.mat,
                          fitted(res.comp)[rm.idx])

knitr::kable(round(valid.pred[,c(1,3,2)], 2))

## ----prediction-exercise-2, purl = TRUE, message = FALSE, echo = 1:4, results = 'asis'----
rm.idx <- sample(nrow(dat), nrow(dat)/10)
dat.cv <- dat
dat.cv[rm.idx, 'phenotype'] <- NA
## Re-fit the model and build table
res.comp.cv <-
  remlf90(fixed   = phenotype ~ 1,
          genetic = list(model = 'competition',
                         pedigree = dat[, 1:3],
                         id = 'self',
                         coord = dat[, c('x', 'y')],
                         competition_decay = 1,
                         pec = list(present = TRUE)),
          spatial = list(model = 'AR', 
                         coord = dat[, c('x', 'y')],
                         rho   = c(.3, .8)),
          data = dat.cv,
          method = 'em')

var.comp <- 
  data.frame(
    'Fully estimated' = with(res.comp$var,
                             round(c(genetic[1, 1], genetic[2, 2],
                                     genetic[1, 2]/sqrt(prod(diag(genetic))),
                                     spatial, pec, Residual), digits = 2)),
    'CV estimated' = with(res.comp.cv$var,
                          round(c(genetic[1, 1], genetic[2, 2],
                                  genetic[1, 2]/sqrt(prod(diag(genetic))),
                                  spatial, pec, Residual), digits = 2)),
    row.names = c('direct', 'compet.', 'correl.',
                  'spatial', 'pec', 'residual')
    )

knitr::kable(var.comp)

## ----prediction-exercise-3, purl = TRUE----------------------------------
true.exp.cv <- with(dat[rm.idx, ], phenotype - resid)
round(sqrt(mean((fitted(res.comp.cv)[rm.idx] - true.exp.cv)^2)), 2)

## ----multitrait-fit------------------------------------------------------
## Filter site and select relevant variables
dat <- 
  droplevels(
    douglas[douglas$site == "s3",
            names(douglas)[!grepl("H0[^4]|AN|BR|site", names(douglas))]]
  )

res <- 
  remlf90(
    fixed = cbind(H04, C13) ~ orig,
    genetic = list(
      model = 'add_animal', 
      pedigree = dat[, 1:3],
      id = 'self'),
    data = dat
  )

## ----multitrait-summary, echo = FALSE------------------------------------
summary(res)

## ----multitrait-genetic-covariances--------------------------------------
res$var[["genetic", "Estimated variances"]]

## Use cov2cor() to compute correlations
cov2cor(res$var[["genetic", "Estimated variances"]])

## ----multitrait-fixef-ranef----------------------------------------------
fixef(res)       ## printed in tabular form, but...
unclass(fixef(res))  ## actually a matrix of estimates with attribute "se"

str(ranef(res))
head(ranef(res)$genetic)

## ----multitrait-blups----------------------------------------------------
head(model.matrix(res)$genetic %*% ranef(res)$genetic)

## ----multitrait-var-ini-spec, eval = FALSE-------------------------------
#  initial_covs <- list(
#    genetic = 1e3*matrix(c(1, .5, .5, 1), nrow = 2),
#    residual = diag(2)   # no residual covariances
#  )
#  res <-
#    remlf90(
#      fixed = cbind(H04, C13) ~ orig,
#      genetic = list(
#        model = 'add_animal',
#        pedigree = dat[, 1:3],
#        id = 'self',
#        var.ini = initial_covs$genetic),
#      data = dat,
#      var.ini = list(residual = initial_covs$residual)
#    )

## ----breedR-options------------------------------------------------------
breedR.getOption()

