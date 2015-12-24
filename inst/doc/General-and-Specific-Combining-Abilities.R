## ----factorial-setup, warning = FALSE------------------------------------
## Setup
library(breedR)
library(ggplot2)
set.seed(123)

## Simulation parameters
n.parents <- c(male = 15L,
               female = 10L)
n.matings <- 100
n.replicates <- 40
mu = 10        # Intercept
sigma2_g <- 6  # Genetic variance of the base population
sigma2_s <- 1  # Variance of the SCA
sigma2_e <- 1  # Residual variance
## Generate all crosses and sample a subset
parents.codes <- list(male = seq.int(n.parents['male']),
                      female = n.parents['male'] + seq.int(n.parents['female']))
matings <- expand.grid(parents.codes)
matings <- matings[sample(prod(n.parents), n.matings),]
rownames(matings) <- with(matings, paste(male, female, sep = 'x'))

## Simulated values
GCA = sapply(do.call('c', parents.codes),
             function(x) rnorm(1, mean = 0, sd = sqrt(sigma2_g)))
SCA = sapply(rownames(matings),
             function(x) rnorm(1, mean = 0, sd = sqrt(sigma2_s)))

## Expected phenotype per family
eta.family <- mu + SCA + (GCA[matings$male] + GCA[matings$female])/2

## Realised Breeding Values in the progeny
## (intra-family variance = half genetic variance)
n.progeny <- n.replicates*n.matings
eta.realised <- eta.family + rnorm(n.progeny, sd = sqrt(sigma2_g/2))

dat <- data.frame(Id = max(sapply(parents.codes, max)) + seq.int(n.progeny),
                  rep = rep(seq.int(n.replicates), each = n.matings),
                  matings,
                  eta.realised,
                  y = eta.realised + rnorm(n.progeny, sd = sqrt(sigma2_e)))

## Define variable for the non-additive SCA
dat <- transform(dat,
                 SCA   = factor(paste(male, female, sep = 'x'),
                                levels = rownames(matings)))

## Printing simulated setting
print(table(dat[, c('male', 'female')]), zero.print = "")
str(dat)

## ----factorial-parents-effects, warning=FALSE----------------------------
## Note that I would like to estimate only **one** GCA effect
## However, currently I need to specify two independent random effects with
## two independent variances, which account in reality for the same thing
res <- remlf90(y ~ 1,
               random = ~ male + female + SCA,
               dat = transform(dat,
                               male   = factor(male),
                               female = factor(female)))

## Here, the effects 'female' and 'male' are both estimating GCA/2
## therefore, their variances are Var(GCA)/4 = sigma_g/4
## So, a point estimator for sigma_g would be:
(sigma_g.est <- 4 * mean(res$var[c('female', 'male'), 1]))
## while the BLUPs
PGCA <- c(ranef(res)$male, ranef(res)$female)

## Check fit
qplot(dat$eta, fitted(res)) + geom_abline(intercept=0, slope=1)
qplot(GCA, PGCA) + geom_abline(intercept=0, slope=1)
qplot(SCA, ranef(res)$SCA) + geom_abline(intercept=0, slope=1)

summary(res)

## ----factorial-pedigree, warning=FALSE-----------------------------------
res.add <- remlf90(y ~ 1,
                   random  = ~ SCA,
                   genetic = list(model    = 'add_animal',
                                  pedigree = dat[, c('Id', 'male', 'female')],
                                  id       = 'Id'),
                   dat = dat)

# Check fit
qplot(dat$eta, fitted(res.add)) + geom_abline(intercept=0, slope=1)
# Predicted GCAs for the parents
PGCA.add <- ranef(res.add)$genetic[do.call('c', parents.codes)]
qplot(GCA, PGCA.add) + geom_abline(intercept=0, slope=1)
# Predicted SCAs for the families
qplot(SCA, ranef(res.add)$SCA) + geom_abline(intercept=0, slope=1)

summary(res)

