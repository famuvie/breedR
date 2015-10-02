old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

### Spatial-Autoregressive model for spatially correlated traits ###

# Eucaliptus Globulus dataset
# Cappa and Cantet (2007)
data(globulus)
cat(paste(comment(globulus), collapse='\n'))

if( !require('ggplot2') ) {
  stop('This demo requires package ggplot2. Please install.')
}

# Variables and structure:
head(globulus)
str(globulus)
  # Notice that 'Factor' variables take a discrete number of values
  # and can be used either as fixed or random effects

# A spatial 'AR' process without specification of the autoregressive parameters
# will automatically fit several models with different combinations.
res.ar  <- remlf90(fixed  = phe_X ~ gg,
                   genetic = list(model = 'add_animal', 
                                  pedigree = globulus[,1:3],
                                  id = 'self'), 
                   spatial = list(model = 'AR', 
                                  coord = globulus[, c('x','y')]), 
                   data = globulus)


# You can visualize the log-likelihood of the models with 
# the tried combinations of autoregressive parameters.
qplot(rho_r, rho_c, fill = loglik, geom = 'tile', data = res.ar$rho)

# Refine the grid around the most likely values
rho.grid <- expand.grid(rho_r = seq(.7, .95, length = 4),
                        rho_c = seq(.7, .95, length = 4))
res.ar  <- remlf90(fixed  = phe_X ~ gg,
                   genetic = list(model = 'add_animal', 
                                  pedigree = globulus[,1:3],
                                  id = 'self'), 
                   spatial = list(model = 'AR', 
                                  coord = globulus[, c('x','y')],
                                  rho = rho.grid), 
                   data = globulus)


# A summary shows the selected model with the most likely combination
# of the autoregressive parameters
summary(res.ar)

### Predicted spatial effect in the observed locations ###
plot(res.ar, 'spatial')

### Prediction also in unobserved locations (full grid) ###
plot(res.ar, 'fullspatial')
