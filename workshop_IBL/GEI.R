## ----libraries, include = FALSE------------------------------------------
library(breedR)
library(tidyr)    # reshaping data: spread() and separate()
library(dplyr)    # pipe operator %>% and data-wrangling: 
                  #   mutate(), group_by(), summarise(), filter()
library(ggplot2)  # implementation of grammar of graphics (Wilkinson L. 2005):
                  # http://www.springer.com/us/book/9780387245447
library(GGally)   # ggpairs()
library(knitr)    # compile Rnw

## ----description-structure} {.smaller------------------------------------
data("douglas")
str(douglas)

## ----description-origxsite-----------------------------------------------
with(douglas, table(orig, site))

## ----description-visualization, fig.width  = 8, warning = FALSE----------
ggplot(douglas, aes(x, y, color = block)) +
  geom_point(show_guide = FALSE) +
  facet_wrap(~ site)

## ----cleanup, message = FALSE--------------------------------------------

## create site-wise block variables (i.e. bl_i = bl * Ind(i))
dat <- transform(douglas, bl1 = block, bl2 = block, bl3 = block)
dat$bl1[dat$site != "s1"] <- NA
dat$bl2[dat$site != "s2"] <- NA
dat$bl3[dat$site != "s3"] <- NA
dat <- droplevels(dat)

## variable family (taken as a maternal effect)
dat$fam <- factor(dat$mum)

## family-site interaction variable
dat <- transform(dat, famxsite = factor(fam:site))

## ----fit-famxsite, message = FALSE---------------------------------------
reml.tree <- remlf90(fixed  = C13 ~ site + orig,
                     random = ~ fam + bl1 + bl2 + bl3 + famxsite,
                     data = dat)

## ----summary, echo = FALSE-----------------------------------------------
summary(reml.tree) 

## ----diagnostics, echo = FALSE, message = FALSE, warning = FALSE---------
data.frame(Residuals = resid(reml.tree)) %>% 
  ggplot(aes(x = Residuals)) +
  geom_histogram()

data.frame(C13 = dat$C13, 
           Prediction = fitted(reml.tree), 
           Site = dat$site) %>% 
  ggplot(aes(x = C13, y = Prediction, color = Site)) + 
  geom_point() +
  geom_abline(int = 0, sl = 1)

## ----plot-interaction-effect, echo = FALSE, fig.width=8------------------
## BLUPs of interaction
dat.fsi <-
  data.frame(
    famxsite = factor(1:nlevels(dat$famxsite),
                      labels = levels(dat$famxsite)),
    pred_fsi = as.vector(ranef(reml.tree)$famxsite)
  ) %>% 
  tidyr::separate(
    famxsite, c('fam', 'site'),
    convert = TRUE,
    remove  = FALSE
  ) %>% 
  transform(
    fam = as.factor(fam),
    site = as.factor(site)
  )

ggplot(dat.fsi, aes(as.numeric(site), pred_fsi, group = fam)) +
  geom_line(color = 'darkgray') +
  scale_x_continuous(breaks = 1:3,
                     labels = paste0('s', 1:3)) + 
  theme(axis.ticks.x = element_blank()) +
  xlab("Site") +
  ylab("Interaction")

## ----select-families-----------------------------------------------------
famxsite.tbl <- table(dat$fam, dat$site) > 0
fam3.idx <- apply(famxsite.tbl, 1, all)
table(fam3.idx)

## ----ecovalence, results = 'hold'----------------------------------------
(dat.ecov <- 
   dat.fsi %>%
   group_by(fam) %>%
   summarise(ssfam =
               sum(pred_fsi**2)) %>%
   filter(fam3.idx) %>% 
   mutate(Ecovalence=ssfam/sum(ssfam)))

## Check
stopifnot(
  all.equal(sum(dat.ecov$Ecovalence), 1)
)

## ----identify-interactive-families, message = FALSE, fig.height=3, fig.show = 'hold', echo = -5----
threshold <- 0.05  # arbitrarily chosen
dat.ecov$hi <- dat.ecov$Ecovalence > threshold
ggplot(dat.ecov, aes(Ecovalence, fill = hi)) +
  geom_histogram(show_guide = FALSE)

cat('Most interactive families: ', with(dat.ecov, as.numeric(fam[hi])))

## ----plot-most-interactive-families, echo = FALSE, fig.width=8-----------
## Plot only families present in all the three sites
plotdat <- dat.fsi[fam3.idx[dat.fsi$fam], ]

ggplot(plotdat,
       aes(as.numeric(site), pred_fsi, group = fam)) +
  geom_line(color = ifelse(plotdat$fam %in% dat.ecov$fam[dat.ecov$hi],
                           'black', 'darkgray')) +
  scale_x_continuous(breaks = 1:3) + 
  theme(axis.ticks.x = element_blank()) +
  xlab("Site") +
  ylab("Interaction")

## ----type-B-formal, message = FALSE--------------------------------------
r_B.fml <- "(G_3_3_1_1-G_7_7_1_1/3)/(G_3_3_1_1+2*G_7_7_1_1/3)"
reml.tree <- remlf90(
  fixed  = C13 ~ site + orig,
  random = ~ fam + bl1 + bl2 + bl3 + famxsite,
  ### <b>
  progsf90.options = paste("se_covar_function r_B", r_B.fml),
  ### </b>
  data = dat
)

## ----type-B-summary, echo = FALSE----------------------------------------
summary(reml.tree)

## ----predicted-genetic-effects, echo = FALSE, fig.width=8----------------
## compute genetic values by family and site
dat.fsi <- transform(
  dat.fsi,
  pred_gen = ranef(reml.tree)$fam[dat.fsi$fam] + pred_fsi
)

## only plot families in all the three sites
plotdat <- dat.fsi[fam3.idx[dat.fsi$fam], ]
ggplot(plotdat, aes(as.numeric(site), pred_gen, group = fam)) +
  geom_line(color = ifelse(plotdat$fam %in% dat.ecov$fam[dat.ecov$hi],
                           'black', 'darkgray')) +
  scale_x_continuous(breaks = 1:3,
                     labels = paste0('s', 1:3)) + 
  theme(axis.ticks.x = element_blank()) +
  xlab("Site") +
  ylab("Interaction")


## ----correlations-reshape, results = 'hold', echo = -12------------------
## reshape the genetic values
## by family

(dat.gen <- 
   
   dat.fsi %>% 
   
   select(fam, site, pred_gen) %>% 
   
   tidyr::spread(site, pred_gen) %>% 
   
   tbl_df()
)


## ----correlations-compute, echo = FALSE, warning = FALSE, fig.width=8----
## plot and pearson correlation
GGally::ggpairs(dat.gen[,-1])

## ----correlations-rank, echo = FALSE-------------------------------------
## Table of kendal and spearman rank correlation measurements
cor.gen <- cor(dat.gen[, -1], method = 'kendall', use = "na.or.complete")
cor.gen[lower.tri(cor.gen)] <-
  cor(dat.gen[, -1],
      method = 'spearman',
      use = "na.or.complete")[lower.tri(cor.gen)]

# print(cor.gen, digits = 3)
kable(cor.gen, digits = 3)

## ----describe-variances, echo = FALSE, warning=FALSE, message=FALSE------
## Empirical variance of phenotype by site
dat %>%
  group_by(fam, site) %>%
  summarise(mean_C13 = mean(C13, na.rm = TRUE)) %>% 
  ggplot(aes(site, mean_C13)) +
  geom_violin()

## Empirical variance of fsi BLUPs by site
# dat.fsi %>% group_by(site) %>% summarise(var_pfsi = var(pred_fsi))

ggplot(dat.fsi, aes(pred_fsi)) + 
  geom_histogram() + 
  facet_grid(site~.) + 
  xlab('Predicted family-site interaction')



## ----cleanup-independent-fsi, message = FALSE, echo = FALSE--------------
## First approximation: three independent effects.

## create site-wise family variables (i.e. fi = fam * Ind(i))
dat <- transform(dat, f1 = fam, f2 = fam, f3 = fam)
dat$f1[dat$site != "s1"] <- NA
dat$f2[dat$site != "s2"] <- NA
dat$f3[dat$site != "s3"] <- NA
# dat <- droplevels(dat)
# non observed values will be estimated as zero

## ----fit-independent-fsi, message = FALSE, echo = FALSE------------------
reml.test <- remlf90(fixed  = C13 ~ site + orig,
                     random = ~ bl1 + bl2 + bl3 + f1 + f2 + f3,
                     data = dat
                     # , breedR.bin = 'submit'
                     )


## ----summary-independent-fsi, echo = FALSE-------------------------------
summary(reml.test)

## ----diagnostics-independent-fsi, echo = FALSE, message = FALSE, warning = FALSE----
ggplot(data.frame(Residuals = resid(reml.test)),
       aes(x = Residuals)) +
  geom_histogram()

ggplot(data.frame(C13 = dat$C13, Prediction = fitted(reml.test)),
       aes(x = C13, y = Prediction, color = dat$site)) + 
  geom_point() +
  geom_abline(int = 0, sl = 1)

## ----plot-PBV-independent-effect, echo = FALSE, fig.width=8--------------
## BLUPs of interaction
# str(dat.fsi)

## recover interaction prediction in the right order
get_pred_fsi <- function(x, fam, site) {
  sapply(seq_along(fam), function(i) x[[site[i]]][fam[i]])
}

dat.fsi$pred_gen_indep <- 
  get_pred_fsi(ranef(reml.test)[paste0('f', 1:3)],
               as.numeric(dat.fsi$fam),
               as.numeric(dat.fsi$site))

ggplot(dat.fsi, aes(as.numeric(site), pred_gen_indep, group = fam)) +
  geom_line(color = ifelse(dat.fsi$fam %in% dat.ecov$fam[dat.ecov$hi],
                           'black', 'darkgray')) +
  scale_x_continuous(breaks = 1:3,
                     labels = paste0('s', 1:3)) + 
  theme(axis.ticks.x = element_blank()) +
  xlab("Site") +
  ylab("Predicted genetic effects")


## ----unified-effect-vs-independent-effects, echo = FALSE, fig.width=6, fig.height=5----
ggplot(dat.fsi, aes(x = pred_gen, y = pred_gen_indep, color = site)) + 
  xlab('Single interaction model') +
  ylab('Independent interaction effects') +
  geom_point() + 
  geom_abline(int = 0, sl = 1)

## ----correlations-reshape-indep, warning=FALSE, echo = FALSE, fig.width=8, fig.height=5----
dat.gen <- 
  dat.fsi %>% 
  select(fam, site, pred_gen_indep) %>% 
  tidyr::spread(site, pred_gen_indep) %>% 
  tbl_df()

GGally::ggpairs(dat.gen[,-1])

## ----model-hack, echo = FALSE--------------------------------------------
## To make the three effects correlated we need some hacking.

## Start with the independent-interactions model
effects <- reml.test$effects

## Group the three independent effects to allow covariation
effects$fsi <- breedR:::effect_group(
  list(f1 = effects$f1$effects[[1]],
       f2 = effects$f2$effects[[1]],
       f3 = effects$f3$effects[[1]]),
  cov.ini = matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3, 3)
)

## Remove previous independent effects
effects$f1 <- effects$f2 <- effects$f3 <- NULL

## Manually build pf90 object
pf90 <- breedR:::progsf90(model.frame(reml.test), 
                          effects, opt = c("sol se"), 
                          res.var.ini = 1)

tmpdir <- tempdir()

breedR:::write.progsf90(pf90, dir = tmpdir)

cdir <- setwd(tmpdir)

reml.out <- system2(file.path(breedR.getOption("breedR.bin"), 'airemlf90'), 
                    input  = 'parameters',
                    stdout = TRUE)

reml.test2 <- breedR:::parse_results(file.path(tmpdir, 'solutions'), effects, model.frame(reml.test), reml.out, 'ai', reml.test$call)
class(reml.test2) <- c('breedR', 'remlf90')

setwd(cdir)

## ----summary-correlated, echo = FALSE------------------------------------
summary(reml.test2)
vm <- reml.test2$var$fsi
sdm <- diag(1/sqrt(diag(vm)))
round(sdm %*% vm %*% sdm, 3)

## ----diagnostics-correlated-fsi, echo = FALSE, message = FALSE, warning = FALSE----
ggplot(data.frame(Residuals = resid(reml.test2)),
       aes(x = Residuals)) +
  geom_histogram()

ggplot(data.frame(C13 = dat$C13, Prediction = fitted(reml.test2)),
       aes(x = C13, y = Prediction, color = dat$site)) + 
  geom_point() +
  geom_abline(int = 0, sl = 1)

## ----plot-PBV-correlated-effect, echo = FALSE, fig.width=8---------------
## recover interaction prediction in the right order
get_pred_fsi <- function(x, fam, site) {
  sapply(seq_along(fam), function(i) x[[site[i]]][fam[i]])
}

dat.fsi$pred_gen_cor <- 
  get_pred_fsi(ranef(reml.test2)[paste0('f', 1:3)],
               as.numeric(dat.fsi$fam),
               as.numeric(dat.fsi$site))

ggplot(dat.fsi, aes(as.numeric(site), pred_gen_cor, group = fam)) +
  geom_line(color = ifelse(dat.fsi$fam %in% dat.ecov$fam[dat.ecov$hi],
                           'black', 'darkgray')) +
  scale_x_continuous(breaks = 1:3,
                     labels = paste0('s', 1:3)) + 
  theme(axis.ticks.x = element_blank()) +
  xlab("Site") +
  ylab("Interaction")

## ----unified-effect-vs-correlated-effects, echo = FALSE, fig.width=6, fig.height=5----
ggplot(dat.fsi, aes(x = pred_gen, y = pred_gen_cor, color = site)) + 
  xlab('Single interaction model') +
  ylab('Correlated interaction effects') +
  geom_point() + 
  geom_abline(int = 0, sl = 1)


## ----correlations-reshape-cor, warning=FALSE, echo = FALSE, fig.width=8, fig.height=5----
dat.gen <- 
   dat.fsi %>% 
   select(fam, site, pred_gen_cor) %>% 
   tidyr::spread(site, pred_gen_cor) %>% 
   tbl_df()

GGally::ggpairs(dat.gen[,-1])

## ----prepare-famxsite-splines--------------------------------------------
## Use breedR internal functions to compute splines models on each site
sp1 <- breedR:::breedr_splines(dat[dat$site == 's1', c('x', 'y')])
sp2 <- breedR:::breedr_splines(dat[dat$site == 's2' & !is.na(dat$x) & !is.na(dat$y), c('x', 'y')])
sp3 <- breedR:::breedr_splines(dat[dat$site == 's3', c('x', 'y')])

## Manually build the full incidence matrices with 0
## and the values computed before
mm1 <- model.matrix(sp1)
inc.sp1 <- Matrix::Matrix(0, nrow = nrow(dat), ncol = ncol(mm1))
inc.sp1[dat$site == 's1', ] <- mm1

mm2 <- model.matrix(sp2)
inc.sp2 <- Matrix::Matrix(0, nrow = nrow(dat), ncol = ncol(mm2))
inc.sp2[dat$site == 's2' & !is.na(dat$x) & !is.na(dat$y), ] <- mm2

mm3 <- model.matrix(sp3)
inc.sp3 <- Matrix::Matrix(0, nrow = nrow(dat), ncol = ncol(mm3))
inc.sp3[dat$site == 's3', ] <- mm3


## ----fit-famxsite-splines, message=FALSE---------------------------------
reml.spl <- remlf90(
  fixed  = C13 ~ site + orig,
  random = ~ fam + famxsite,
  generic = list(sp1 = list(inc.sp1,
                            breedR:::get_structure(sp1)),
                 sp2 = list(inc.sp2,
                            breedR:::get_structure(sp2)),
                 sp3 = list(inc.sp3,
                            breedR:::get_structure(sp3))),
  data = dat,
  method = 'em'
)

## ----plot-splines, fig.show='hold', echo = FALSE-------------------------
breedR:::spatial.plot(
  dat = data.frame(
    coordinates(sp1), 
    z = as.vector(mm1 %*% ranef(reml.spl)$sp1)
  )
)
breedR:::spatial.plot(
  dat = data.frame(
    coordinates(sp2), 
    z = as.vector(mm2 %*% ranef(reml.spl)$sp2)
  )
)
breedR:::spatial.plot(
  dat = data.frame(
    coordinates(sp3), 
    z = as.vector(mm3 %*% ranef(reml.spl)$sp3)
  )
)

