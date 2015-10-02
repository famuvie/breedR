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
  geom_abline(int = 0,
              sl = 1,
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
  geom_abline(int = 0,
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

## ----spatial-exercise----------------------------------------------------
rho.grid <- expand.grid(rho_r = seq(.7, .95, length = 4),
                        rho_c = seq(.7, .95, length = 4))

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
  geom_abline(int = 0, sl = 1, col = 'darkgray') + 
  geom_text(label = labels, parse = TRUE, hjust = -.5)


## ----generic-example, message = FALSE------------------------------------
## Fit a blocks effect using generic
inc.mat <- model.matrix(~ 0 + bl, globulus)
cov.mat <- diag(nlevels(globulus$bl))
res.blg <- remlf90(fixed  = phe_X ~ gg,
                   generic = list(block = list(inc.mat,
                                               cov.mat)),
                   data   = globulus)

## ----summary-generic-----------------------------------------------------
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
pred.BV.loo.mat <- with(ranef(res.comp.loo), cbind(`genetic_direct`,
                                               `genetic_competition`))
pred.genetic.loo <- Zd[rm.idx, ] %*% pred.BV.loo.mat

valid.pred <- 
  data.frame(True = with(dat[rm.idx, ],
                         c(BV1, BV2, rm.exp)),
             Pred.loo = c(pred.genetic.loo,
                          fitted(res.comp.loo)[rm.idx]),
             row.names = c('direct BV', 'competition BV', 'exp. phenotype'))

knitr::kable(round(valid.pred, 2))

## ----breedR-options------------------------------------------------------
breedR.getOption()

