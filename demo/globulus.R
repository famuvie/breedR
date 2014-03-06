## Eucaliptus Globulus dataset
## Cappa and Cantet (2007)

data(globulus)
cat(paste(comment(globulus), collapse='\n'))

# Variables and structure:
head(globulus)
str(globulus)
  # Notice that 'Factor' variables take a discrete number of values
  # and can be used either as fixed or random effects

### Blocks ###
# Fit a model with genetic group as a fixed effect, block as a random effect
# and the pedigree-based additive genetic effect.
res.blk <- remlf90(fixed  = phe_X ~ gg,
                     random = ~ bl,
                     genetic = list(model = 'add_animal', 
                                    var.ini = 10, 
                                    pedigree = globulus[,1:3],
                                    id = 'self'), 
                     data = globulus,
                     method = 'em')

### Splines ###
# Alternativelly, account for the environmental variation with a
# continuous spatial effect instead of using discrete blocks.
res.spl  <- remlf90(fixed  = phe_X ~ gg,
                     genetic = list(model = 'add_animal', 
                                    var.ini = 10, 
                                    pedigree = globulus[,1:3],
                                    id = 'self'), 
                     spatial = list(model = 'Cappa07', 
                                    coord = globulus[, c('x','y')], 
                                    knots = c(7, 7),
                                    var.ini = 30), 
                     data = globulus,
                     method = 'em')

### AR1 x AR1 ###
# A further spatial approach with a separable First order Autoregressive
# process on the rows and colums.
# You can fix the values of the autocorrelation parameters, 
# or let the program try several combinations and select the most likely.
res.ar  <- remlf90(fixed  = phe_X ~ gg,
                   genetic = list(model = 'add_animal', 
                                  var.ini = 10, 
                                  pedigree = globulus[,1:3],
                                  id = 'self'), 
                   spatial = list(model = 'AR', 
                                  coord = globulus[, c('x','y')],
                                  rho = c(.85, .8), 
                                  var.ini = 30), 
                   data = globulus,
                   method = 'em')

# Summaries
summary(res.blk)
summary(res.spl)
summary(res.ar)
  # Notice how the continuous spatial model takes longer but gathers more
  # variability by reducing the residual and genetic variances a bit.



### Comparison of spatial effects ###

# Spatial effects estimates in the positions of the observations
# Ordered by blocks and increasing values of the splines effect
spatial.dat <- transform(globulus,
                         Block = ranef(res.blk)$bl$value[bl],
                         Splines = res.spl$spatial$fit$z,
                         AR1xAR1 = res.ar$spatial$fit$z)
ord <- with(spatial.dat, order(Block, Splines))

ggplot(cbind(melt(spatial.dat, id = 1:9), Ind = order(ord)),
       aes(Ind, value)) + 
  geom_point(aes(col = variable))


# Plot the Block effects
(p.bl <- ggplot(transform(spatial.dat,
                         z = Block,
                         model = "Block"),
               aes(x, y)) +
  geom_tile(aes(fill = z)) +
  coord_fixed() +
  facet_wrap(~ model) +
  scale_fill_gradient(low='green', high='red'))

# For models with specific spatial effects, breedR provides convenience methods.
plot(res.spl)
plot(res.ar)


# We can compare the results under the same scale with the convenience function 
# compare.plots()
compare.plots(list(Block = p.bl,
                   Splines = plot(res.spl),
                   AR1xAR1 = plot(res.ar)))



### Prediction in unobserved locations ###

# Unfortunately, we can't compare.plots() two tiles with different resolutions.
# We need to manually set up a common scale and grid.arrange() the plots
scale.limits <- range(c(res.spl$spatial$prediction$z, res.ar$spatial$prediction$z))
pred.plots <- list(
  qplot(x, y, fill=z, geom='tile', data = transform(res.spl$spatial$prediction,
                                                    model = 'Splines')) + 
    scale_fill_gradient(low='green', high='red', limits = scale.limits) +
    coord_fixed() + theme(legend.position = "none") + facet_wrap( ~ model),
  qplot(x, y, fill=z, geom='tile', data = transform(res.ar$spatial$prediction,
                                                    model = 'AR1xAR1')) + 
    scale_fill_gradient(low='green', high='red', limits = scale.limits) +
    coord_fixed()  + facet_wrap( ~ model))
require(gridExtra)
grid.arrange(pred.plots[[1]], pred.plots[[2]], ncol = 2)
