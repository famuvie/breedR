## Eucaliptus Globulus dataset
## Cappa and Cantet (2007)

data(globulus)
cat(paste(comment(globulus), collapse='\n'))

# Variables and structure:
head(globulus)
str(globulus)
  # Notice that 'Factor' variables take a discrete number of values
  # and can be used either as fixed or random effects

# Fit a model with genetic group as a fixed effect, block as a random effect
# and the genetic effect related with the pedigree.
res.block <- remlf90(fixed  = phe_X ~ gg,
                     random = ~ bl,
                     genetic = list(model = 'add_animal', 
                                    var.ini = 10, 
                                    pedigree = globulus[,1:3],
                                    id = 'self'), 
                     data = globulus,
                     method = 'em')

# Alternativelly, account for the environmental variation with a
# continuous spatial effect instead of using discrete blocks
res.spat  <- remlf90(fixed  = phe_X ~ gg,
                     genetic = list(model = 'add_animal', 
                                    var.ini = 10, 
                                    pedigree = globulus[,1:3],
                                    id = 'self'), 
                     spatial = list(model = 'Cappa07', 
                                    coord = globulus[, c('x','y')], 
                                    var.ini = 30), 
                     data = globulus,
                     method = 'em')

# Summaries
summary(res.block)
summary(res.spat)
  # Notice how the continuous spatial model takes a longer but gathers more
  # variability by reducing the residual and genetic variances a bit.


# Random effect estimates in the positions of the observations
str(ranef(res.spat))

# Plotting spatial effect
# (not only in the observation locations, but in the whole region)
qplot(x, y, fill=z, geom='tile', data = res.spat$spatial$prediction) + 
  scale_fill_gradient(low='green', high='red')

