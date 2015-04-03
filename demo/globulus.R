old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

## Eucaliptus Globulus dataset
## Cappa and Cantet (2007)
if( !require('ggplot2') ) {
  stop('This demo requires package ggplot2. Please install.')
}
data(globulus)
cat(paste(comment(globulus), collapse='\n'))

# Variables and structure:
head(globulus)
str(globulus)
  # Notice that 'Factor' variables take a discrete number of values
  # and can be used either as fixed or random effects

# The genetic component
gen.globulus <- list(model    = 'add_animal', 
                     pedigree = globulus[,1:3],
                     id       = 'self')

### Blocks ###
# Fit a model with genetic group as a fixed effect, block as a spatial random
# effect and the pedigree-based additive genetic effect.
res.blk <- remlf90(fixed  = phe_X ~ gg,
                   genetic = gen.globulus, 
                   spatial = list(model = 'blocks', 
                                  coord = globulus[, c('x','y')],
                                  id = 'bl'),
                   data = globulus)

### Splines ###
# Alternativelly, account for the environmental variation with a
# continuous spatial effect instead of using discrete blocks.
# We use 'em' method as AI-REML doesn't mix well with splines.
res.spl  <- remlf90(fixed  = phe_X ~ gg,
                    genetic = gen.globulus, 
                    spatial = list(model = 'splines', 
                                   coord = globulus[, c('x','y')], 
                                   n.knots = c(7, 7)), 
                    data = globulus,
                    method = 'em')

### AR1 x AR1 ###
# A further spatial approach with a separable First order Autoregressive
# process on the rows and colums.
# You can fix the values of the autocorrelation parameters, 
# or let the program try several combinations and select the most likely.
res.ar  <- remlf90(fixed  = phe_X ~ gg,
                   genetic = gen.globulus, 
                   spatial = list(model = 'AR', 
                                  coord = globulus[, c('x','y')],
                                  rho = c(.85, .8)), 
                   data = globulus)

# Summaries
summary(res.blk)
summary(res.spl)
summary(res.ar)
  # Notice how the continuous spatial model takes longer but gathers more
  # variability by reducing the residual and genetic variances a bit.


### Comparison of spatial effects ###

# Spatial effects estimates in the positions of the observations
# multiply the spatial incidence matrix by the BLUP of the spatial random effect
space.pred <- function(x) model.matrix(x)$random$spatial %*% ranef(x)$spatial

# Ordered by blocks and increasing values of the splines effect
spatial.dat <- transform(globulus,
                         Blocks  = space.pred(res.blk),
                         Splines = space.pred(res.spl),
                         AR1xAR1 = space.pred(res.ar))
ord <- with(spatial.dat, order(Blocks, Splines))

ggplot(cbind(reshape2::melt(spatial.dat, id = 1:9), Ind = order(ord)),
       aes(Ind, value)) + 
  geom_point(aes(col = variable))


# For models with specific spatial effects, breedR provides convenience plotting
# methods.
(p.blk <- plot(res.blk, type = 'spatial'))
(p.spl <- plot(res.spl, type = 'spatial'))
(p.ar  <- plot(res.ar,  type = 'spatial'))


# We can compare the results under the same scale with the convenience function 
# compare.plots()
compare.plots(list(Blocks  = p.blk,
                   Splines = p.spl,
                   AR1xAR1 = p.ar))


### Prediction of the spatial effect in unobserved locations ###

# We can use the 'fullspatial' plot type instead
compare.plots(list(Blocks  = plot(res.blk, type = 'fullspatial'),
                   Splines = plot(res.spl, type = 'fullspatial'),
                   AR1xAR1 = plot(res.ar, type = 'fullspatial')))



