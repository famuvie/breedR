## Plotting and diagnostics tools in breedR

data(globulus)

res  <- remlf90(fixed  = phe_X ~ gg,
                genetic = list(model = 'add_animal', 
                               pedigree = globulus[,1:3],
                               id = 'self'), 
                spatial = list(model = 'AR', 
                               coord = globulus[, c('x','y')],
                               rho = c(.8, .9)), 
                data = globulus)

# Visualize the phenotype, model fit, spatial component or residuals, 
(p.phe <- plot(res, type = 'phenotype'))
(p.fit <- plot(res, type = 'fitted'))
plot(res, type = 'spatial')
plot(res, type = 'residuals')

# Compare observed and fitted values side by side in the same scale
compare.plots(list(Phenotype = p.phe,
                   Fitted    = p.fit))

# Compare the spatial and residual components
compare.plots(list(Spatial = plot(res, type = 'sp'),
                   Residual= plot(res, type = 're')))

# (semi)Variograms of the residuals
# By default, isotropic, row/column anisotropy (3d surface and heat map)
# and fully anisotropic
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
