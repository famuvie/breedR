# We use a simulated dataset with known parameters
# which we try to infer from phenotype observations
data(m1)
# It has been simulated with the Metagene program.
# The package has specific methods for 'metagene' objects
summary(m1)
plot(m1)

# Fit Mixed Model using EM-REML
# The formula specifies fixed effects 
# (random unstructured effects are not yet implemented)
# the genetic and spatial effects are optional
# ~~~ This takes a few minutes !!! ~~~

dat <- as.data.frame(m1)

res.f90 <- remlf90(formula = phe_X ~ sex, 
                   genetic = list(model = 'add_animal', 
                                  var.ini = 10, 
                                  pedigree = get_pedigree(m1),
                                  id = 'self'), 
                   spatial = list(model = 'Cappa07', 
                                  coord = coordinates(m1), 
                                  var.ini = 300), 
                   data = dat,
                   method = 'em')

#### Summary of results ####
summary(res.f90)

FO  <- fitted(res.f90)          # Fitted Observations
PBV <- ranef(res.f90)$genetic   # Predicted Breeding Values
PSE <- ranef(res.f90)$spatial   # Predicted Spatial Effect

# Fitted values vs. Observed phenotypes by generation
qplot(phe_X, FO, color = sex, data = dat) +
  geom_abline(int=0, slope=1)

# Predicted genetic values vs. true Breeding values by generation
qplot(BV_X, PBV, color = sex, data = dat) +
  geom_abline(int = 0, slope = 1, col = 'gray')

# Linear correlation between true and fitted breeding values
cor(dat$BV_X, PBV)  # 0.90 !!


#### Spatial structure ####

# True generating surface
plot(m1, 'spatial')

# Estimated surface 
# (not only in the observation locations, but in the whole region)
qplot(irow, icol, fill = PSE, geom = 'tile', data = dat) + 
  scale_fill_gradient(low = 'green', high = 'red')

# The splines interpolation smoothes out many details.
# Globally, however, the predictions are around the true values:

# Predicted Spatial effect vs. true spatial structure
qplot(sp_X, PSE, color = icol, data = dat) +
  geom_abline(int = 0, slope = 1, col = 'gray')

# Do the differences have any spatial pattern?
# They shouldn't:
qplot(irow, icol, fill = sp_X - PSE, geom = 'tile', data = dat) + 
  scale_fill_gradient(low = 'green', high = 'red')
