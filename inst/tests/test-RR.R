### Random Regression testing pilot ###

## Simulate some RR data (following Shaeffer's 2002 notation)
## y_ikn:t = mu + F_i + r(a, t, m)_k + e_ikn:t
## where 
##   F_i: fixed effect independent of the longitudinal scale
##   r(a, x, m)_k = \sum_{l=0}^m a_kl x_ik:l is the regression function
##        for the BV of the animal k. This is a polynomial of order m.
##        It depends on (m+1) independent coefficients a_kl.

                                    # Numbers of :
N <- list(fnd = c(10, 10),          # founders
          ind = 1e3,                # individuals (descendants)
          sam = 5:10,               # number of samples per individual (range)
          lev = 3)                  # levels of the fixed effect

                                    # Variances for
sigma2 <- list(a = 2,               # additive direct genetic effect
               e = 1)               # residual

ord      <- 4                  # polynomial order for the RR function
F_values <- 2^(1:N$lev - 1)    # Values of the levels of the fixed effect
t_range  <- c(0, 10)           # Range of the longitudinal variable

# Random pedigree (random mating of founders)
ped <- breedR.sample.pedigree(N$ind, N$fnd, check.factorial = FALSE)

# Simulate independent BV coefficients (including founders)
# coefficients are independent for the same individual, but are
# corelated between individuals
true_BV <- breedR.sample.BV(ped, sigma2$a, N = ord + 1)

# Times of observation for each individual
# a random number of times (N$sam) uniformly sampled from the specified t_range
tobs <- replicate(N$ind, sort(runif(sample(N$sam, 1), t_range[1], t_range[2])))

# Contribution of additive genetic effect to phenotype
ind_obs <- rep(sum(N$fnd) + 1:N$ind, times = sapply(tobs, length))
bv_obs <- true_BV[ind_obs,]
bv_coef <- legendre_values(unlist(tobs), ord)
bv_cont <- rowSums(bv_obs * bv_coef)

# Simulated data.frame
dat <- transform(data.frame(ind = ind_obs,
                            F   = factor(sample(LETTERS[1:N$lev],
                                                length(ind_obs),
                                                replace = TRUE)),
                            tobs = unlist(tobs)),
                 phe = F_values[as.numeric(F)] + bv_cont + rnorm(length(ind)))


########## End of simulation ##########





########## Pilot demonstration of model fit with breedR ############

# Input: pedigree (ped) and observations data.frame (dat)

# Start by fitting the model as if we had one observation per individual
# (i.e. as if it wasn't longitudinal)
# This is only in order to get a first sketch of the model structure

res_dummy <- remlf90(phe ~ 1 + F,
                     genetic = list(model = 'add_animal',
                                    pedigree = ped,
                                    id = 'ind'),
                     data = dat)


# Now we are going to "manually" modify the relevant internal structures

# The 'genetic' effect will no longer be one value, but (ord + 1)
new_eff <- within(res_dummy$effects, {
  # positions in data of poly values
  genetic$pos <- genetic$pos + 0:ord
  # levels of each regressor
  genetic$levels <- rep(genetic$levels, ord + 1)
  # ord + 1 covariates nested in animal
  genetic$type <- rep(paste('cov',
                            # new position for the animal id
                            res_dummy$effects$genetic$pos+ord+1),
                      ord + 1)
  # covariance matrix for ord + 1 new effects
  genetic$var <- diag(genetic$var, ord + 1)
})

# Change incidence matrix of genetic effec
# Use new breedR function to compute the incidence matrix
inc_mat <- legendre_values(dat$tobs, ord)
colnames(inc_mat) <- paste('a', 1:(ord+1), sep = '_')
new_eff$genetic$gen$B <- cbind(inc_mat, new_eff$genetic$gen$B)

# Build the progsf90 object to be written and passed to Misztal
# Note the residual variance is as in the default. Change if necessary.
pf90 <- breedR:::progsf90(res_dummy$mf, new_eff, opt = c("sol se"), res.var.ini = 1)

# Modify files to be passed to Misztal
breedR:::write.progsf90(pf90, dir = tempdir())

# Run remlf90 (change the script for airemlf90)
binary.path <- breedR:::breedR.call.builtin()
cdir <- setwd(tempdir())
reml.out <- system2(file.path(binary.path, 'remlf90'), 
                    input  = 'parameters',
                    stdout = TRUE)
# Parse solutions
res <- breedR:::parse_results('solutions', new_eff, res_dummy$mf, reml.out, 'em', res_dummy$call)
class(res) <- c('breedR', 'remlf90') 
setwd(cdir)

####################### End of pilot fitting ########################





##################### Diagnostics ########################

summary(res)

# predicted coefficients of breeding values
pred_BV <- sapply(ranef(res)$genetic, function(x) x$value)

# true vs. predicted coefficients
ggplot(data.frame(gen = factor(rep(c('founder', 'progeny'),
                                   times = c(sum(N$fnd), N$ind))),
                  coef = factor(rep(1:(ord+1), each = sum(N$fnd) +N$ind)),
                  true_BV = as.vector(true_BV),
                  pred_BV = as.vector(pred_BV)),
       aes(true_BV, pred_BV)) + 
  geom_point(aes(col = coef)) + 
  geom_abline(int = 0, sl = 1, col = 'darkgray') +
  coord_fixed()



# fitted vs. observations
ggplot(data.frame(fit = fitted(res),
                  obs = dat$phe),
       aes(fit, obs)) + 
  geom_point() + 
  geom_abline(int = 0, sl = 1, col = 'darkgray') + 
  coord_fixed()


# variance components
ggplot(data.frame(name = c(paste('a', 1:(ord+1), sep ='_'), 'resid'),
                  est = with(res$var, c(diag(genetic), res=Residual)),
                  true = c(rep(sigma2$a, ord+1), sigma2$e)),
       aes(true, est)) + 
  geom_point(aes(col = name)) + 
  geom_abline(int = 0, sl = 1, col = 'darkgray') +
  coord_fixed()


# Some BVs in time
n_show = 3
idx <- sample(sum(N$fnd)+1:N$ind, n_show)
x <- seq(t_range[1], t_range[2], length = 101)
true_BVf <- legendre_values(x, ord) %*% t(true_BV[idx,])
pred_BVf <- legendre_values(x, ord) %*% t(pred_BV[idx,])
BVf <- transform(x,
                 t = x,
                 ind = factor(rep(idx, each = length(x), times = 2)),
                 type = factor(rep(c('true', 'pred'), each = length(x)*n_show)),
                 BV = c(as.vector(true_BVf), as.vector(pred_BVf)))

ggplot(BVf, aes(t, BV)) + 
  geom_line(aes(col = ind, lty = type))



#####################################################################
# Conclusions:
# 1. Seems to work fine. It fits the data very well and predicts
#    the BV coefficients and curves nicely
# 2. However, the variance components of the genetic effects are almost
#    doubled!! I'm puzzled with this, since everything else fits well.


#######################################################################
# Notes:
# 1. If in res_dummy$effects are further effects after the 'genetic'
#    you need to shift their position(s) 'pos' by ord

