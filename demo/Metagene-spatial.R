# We use a simulated dataset with known parameters
# which we try to infer from phenotype observations
data(m1)
# It has been simulated with the Metagene program.
# The package has specific methods for 'metagene' objects
summary(m1)
plot(m1)

if( !require('ggplot2') ) {
  stop('This demo requires package ggplot2. Please install.')
}

# Fit Mixed Model using EM-REML
# The formula specifies fixed effects 
# the genetic and spatial effects are optional
# ~~~ This takes a few minutes !!! ~~~

dat <- as.data.frame(m1)

res.f90 <- remlf90(fixed   = phe_X ~ sex, 
                   genetic = list(model = 'add_animal', 
                                  pedigree = get_pedigree(m1),
                                  id = 'self'), 
                   spatial = list(model = 'splines', 
                                  coord = coordinates(m1)), 
                   data = dat,
                   method = 'em')

#### Summary of results ####
summary(res.f90)

FO  <- fitted(res.f90)          # Fitted Observations
PBV <- res.f90$genetic$fit      # Predicted Breeding Values
PSE <- res.f90$spatial$fit$z    # Predicted Spatial Effect

# Fitted values vs. Observed phenotypes by sex
qplot(phe_X, FO, color = sex, data = dat) +
  geom_abline(int=0, slope=1)

# Predicted genetic values vs. true Breeding values by sex
qplot(BV_X, PBV, color = sex, data = dat) +
  geom_abline(int = 0, slope = 1, col = 'gray')

# Linear correlation between true and fitted breeding values
cor(dat$BV_X, PBV)  # 0.84 !!


#### Spatial structure ####

# True generating surface
plot(m1, 'spatial')

# Estimated surface 
# The splines interpolation smoothes out many details.
plot(res.f90, 'spatial')

# Globally, however, the predictions are around the true values:
# Predicted Spatial effect vs. true spatial structure
qplot(sp_X, PSE, data = dat) +
  geom_abline(int = 0, slope = 1, col = 'gray')

# Do the differences have any spatial pattern?
# They shouldn't:
plot(res.f90, z = dat$sp_X - PSE)

# Notice that there is some shorter-range spatial pattern
# It might be worth trying with more knots, but be careful!
# it increases computation time dramatically!
# From the summary, you can see that the previous model
# was fitted with 6 x 6 knots.
# Simply add something like n.knots = c(8, 8) to the list 
# of spatial parameters.
