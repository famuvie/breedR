## ----libraries, include = FALSE------------------------------------------
library(breedR)
library(tidyr)    # reshaping dat: spread() and separate()
library(dplyr)    # pipe operator %>% and dat-wrangling: 
                  #   mutate(), group_by(), summarise(), filter()
library(ggplot2)  # implementation of grammar of graphics (Wilkinson L. 2005):
                  # http://www.springer.com/us/book/9780387245447
library(GGally)   # ggpairs()
library(knitr)    # compile Rnw
library(viridis)

## ----description-structure} {.smaller------------------------------------
data("douglas")
str(douglas)
dat <- droplevels(subset(douglas, site == "s3"))

## ----phenotype, message = FALSE, fig.show='hold'-------------------------
ggplot(dat, aes(C13)) + geom_histogram()

ggplot(dat, aes(x, y, fill = C13)) +
  geom_raster(show_guide = FALSE) +
  coord_fixed() + 
  scale_fill_viridis()

## ----vgram-phenotype, message=FALSE, fig.width=8, fig.height=5-----------
variogram(coord = dat[, c('x','y')], z = dat$C13)

## ----fam_block-fit, message=FALSE----------------------------------------
dat$fam <- factor(dat$mum)
h2.fml <- "4*G_3_3_1_1/(G_2_2_1_1+G_3_3_1_1+R_1_1)"

res.base <- remlf90(
  fixed = C13 ~ orig,
  random = ~ block + fam,
  progsf90.options = paste("se_covar_function h2", h2.fml),
  dat = dat)

## ----fam_block-summary, echo = FALSE-------------------------------------
summary(res.base)

## ----ind_block-fit, message=FALSE, warning=FALSE-------------------------
res.blk <- remlf90(
  fixed = C13 ~ orig,
  genetic = list(model = 'add_animal',
                 pedigree = dat[, c('self','dad','mum')],
                 id = 'self'),
  spatial = list(model = 'blocks',
                 coord = dat[, c('x','y')],
                 id = "block"),
  dat = dat)


## ----ind_block-summary, echo = FALSE-------------------------------------
summary(res.blk)

## ----splines-fit, message=FALSE, warning=FALSE---------------------------
res.spl <- remlf90(
  fixed = C13 ~ orig,
  genetic = list(model = 'add_animal',
                 pedigree = dat[, c('self','dad','mum')],
                 id = 'self'),
  spatial = list(model = 'splines',
                 coord = dat[, c('x','y')]),
  dat = dat,
  method = 'em')


## ----splines-summary, echo = FALSE---------------------------------------
summary(res.spl)

## ----AR-fit, message=FALSE, warning=FALSE--------------------------------
res.ar1 <- remlf90(
  fixed = C13 ~ orig,
  # random = ~ block,
  genetic = list(model = 'add_animal',
                 pedigree = dat[, c('self','dad','mum')],
                 id = 'self'),
  spatial = list(model = 'AR',
                 coord = dat[, c('x','y')],
                 rho = c(.8,.8)),
  dat = dat)


## ----AR-summary, echo = FALSE--------------------------------------------
summary(res.ar1)

## ----residuals-comparison, fig.width=8-----------------------------------
coordinates(res.base) <- dat[, c('x', 'y')]
compare.plots(    # preserve scale
  list(`Family/blocks`      = plot(res.base, 'residuals'),
       `Individual/blocks`  = plot(res.blk, 'residuals'),
       `Individual/splines` = plot(res.spl, 'residuals'),
       `Individual/AR1`     = plot(res.ar1, 'residuals'))
)


## ----spatial-comparison, fig.width = 8-----------------------------------
compare.plots(
  list(Blocks = plot(res.blk, type = 'spatial'),
       Splines = plot(res.spl, type = 'spatial'),
       AR1xAR1 = plot(res.ar1, type = 'spatial'))
)

## ----comp_blk-fit, message = FALSE, warning = FALSE----------------------
res.comp <- remlf90(
  fixed = C13 ~ orig,
  genetic = list(
    model = c('comp'),
    pedigree = dat[, c('self','dad','mum')],
    id = 'self',
    coord = dat[, c('x', 'y')],
    competition_decay = 2, # IC decay 1/distance
    pec = TRUE), # envirmonetal compettion effect
  spatial = list(
    model = 'blocks',
    coord = dat[, c('x','y')],
    id = "block"),
  data = dat,
  method = 'em'
)


## ----comp_blk-summary----------------------------------------------------
summary(res.comp)


## ----genetic-correlation-------------------------------------------------
(S <- res.comp$var$genetic)
SD <- sqrt(diag(1/diag(S)))
(SD %*% S %*% SD)[1,2]

## ----comp_AR-fit, message = FALSE, warning = FALSE-----------------------

res.comp.ar1 <- remlf90(
  fixed = C13 ~ orig,
  genetic = list(
    model = c('comp'),
    pedigree = dat[, c('self','dad','mum')],
    id = 'self',
    coord = dat[, c('x', 'y')],
    competition_decay = 2, # IC decay 1/distance
    pec = list(present = TRUE)), #envirmonetal compettion effect
  spatial = list(
    model = 'AR',
    coord = dat[, c('x','y')],
    rho = c(.8,.8)),
  data = dat,
  method = 'em'
)


## ----comp_AR-summary-----------------------------------------------------
summary(res.comp.ar1)


## ----comp_splines-fit, message = FALSE, warning = FALSE------------------

res.comp.spl <- remlf90(
  fixed = C13 ~ orig,
  genetic = list(
    model = c('comp'),
    pedigree = dat[, c('self','dad','mum')],
    id = 'self',
    coord = dat[, c('x', 'y')],
    competition_decay = 2, # IC decay 1/distance
    pec = list(present = TRUE)), #envirmonetal compettion effect
  spatial = list(
    model = 'splines',
    coord = dat[, c('x','y')]),
  data = dat,
  method = 'em'
)


## ----comp_splines-summary------------------------------------------------
summary(res.comp.spl)

## ----model-comparison-table----------------------------------------------

ml <- list(`no comp.` = res.blk, `comp+blk` = res.comp, 
           `comp+ar1` = res.comp.ar1, `comp+spl` = res.comp.spl)
kable(cbind(loglik = sapply(ml, logLik),
            AIC = sapply(ml, function(x) x$fit$AIC)), digits=2)


