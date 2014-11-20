old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

## Plotting and diagnostics tools in breedR

data(globulus)

# Fit a model without accounting for spatial autocorrelation
res0  <- remlf90(fixed  = phe_X ~ gg,
                genetic = list(model = 'add_animal', 
                               pedigree = globulus[,1:3],
                               id = 'self'), 
                data = globulus)

# To plot spatially arranged values, we need to provide the coordinates of the
# observations
coordinates(res0) <- globulus[, c('x','y')]

# Visualize the phenotype and model fit
(p.phe <- plot(res0, type = 'phenotype'))
(p.fit0 <- plot(res0, type = 'fitted'))

# Compare observed and fitted values side by side *in the same scale*
compare.plots(list(Phenotype = p.phe,
                   Fitted    = p.fit0))

  # The model explains very little!

# The residuals looks very autocorrelated
plot(res0, type = 'residuals')

# (semi)Variograms of the residuals
# By default, isotropic, row/column anisotropy (3d surface and heat map)
# and fully anisotropic
variogram(res0)


# There is obviously the need to account for spatial autocorrelation
res  <- remlf90(fixed  = phe_X ~ gg,
                genetic = list(model = 'add_animal', 
                               pedigree = globulus[,1:3],
                               id = 'self'), 
                spatial = list(model = 'AR', 
                               coord = globulus[, c('x','y')],
                               rho = c(.8, .9)), 
                data = globulus)



# Compare observed and fitted values side by side in the same scale
compare.plots(list(Phenotype      = p.phe,
                   Fitted0        = p.fit0,
                   Fitted_spatial = plot(res, type = 'fitted')))

# Compare the spatial and residual components
compare.plots(list(Spatial = plot(res, type = 'sp'),
                   Residual= plot(res, type = 're')))

# Now the variograms are mostly "flat"
variogram(res)

# Individual variograms can be plotted using
# plot = c('all', 'isotropic', 'anisotropic', 'perspective', 'heat', 'none')
variogram(res, plot = 'heat')

# In all cases, the variogram values can be recovered
# and exploited at user's will
vgm <- variogram(res)
str(vgm, 1)
with(vgm$isotropic, plot(distance, variogram, type = 'l'))

# The colours used for spatial plots are customizable
# See ?breedR.setOption
# For example, you can set colours for B&W output
op <- breedR.setOption(col.seq = c('black', 'white'))
plot(res, type = 'fitted')

# Back to default values
breedR.setOption(op)
