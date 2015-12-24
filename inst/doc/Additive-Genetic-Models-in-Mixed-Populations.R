## ----diallel-setup, warning = FALSE--------------------------------------
## Setup
library(breedR)
library(ggplot2)
set.seed(123)

## Simulation parameters
n.founders <- c(E = 9, J = 9)
sigma2 <- c(E = 3, J = 2, resid = 1)

founders <- data.frame(id = c(paste0('E', 1:n.founders['E']),
                              paste0('J', 1:n.founders['J'])),
                       pop.idx = c(rep(1, n.founders['E']),
                                   rep(2, n.founders['J'])),
                       BV = c(rnorm(n.founders['E'], sd = sqrt(sigma2['E'])),
                              rnorm(n.founders['J'], sd = sqrt(sigma2['J']))))


n.obs <- sum(n.founders)*150
obs.parents.idx <- matrix(sample(nrow(founders), 2*n.obs, replace = TRUE), ncol = 2)

## The 'family' is independent of the order of the parents
## While the population only takes into account the origin
## e.g. cross2fam(c('J3', 'E1')) gives 'E1:J3'
## while cross2pop(c('J3', 'E1')) gives 'EJ'
cross2fam <- function(x) paste(founders$id[sort(x)], collapse = ':')
cross2pop <- function(x) paste(names(n.founders)[founders$pop.idx[sort(x)]], collapse = '')

## Mendelian sampling term
msp <- function(x) {
  ss <- sigma2[founders$pop.idx[x]]
  s2 <- (ss[1] + ss[2])/4
  rnorm(1, sd = sqrt(s2))
}

dat <- data.frame(
  id  = sum(n.founders) + seq.int(n.obs),
  dad = obs.parents.idx[, 1],
  mum = obs.parents.idx[, 2],
  fam = apply(obs.parents.idx, 1, cross2fam),
  sp  = apply(obs.parents.idx, 1, cross2pop),
  bv  = apply(obs.parents.idx, 1, 
              function(x) mean(founders$BV[x]) + msp(x)),
  resid = rnorm(n.obs, sd = sqrt(sigma2['resid'])))

dat <- transform(dat,
                 y = bv + resid)

## Printing simulated setting
print(table(dat[, c('mum', 'dad')]), zero.print = "")
str(dat)

## ----overall-genetic-structure-------------------------------------------

## Build a pedigree for the whole mixed population
## and get the kinship matrix A
ped <- build_pedigree(1:3, data = dat)
A <- pedigreemm::getA(ped)

## Build the full incidence matrix
Z <- as(dat$id, 'indMatrix')

## Give the index vector of additive-genetic random effects
## that belong to one subpopulation;
## 'E', 'J' (founders) or 'EE', 'EJ' or 'JJ' (offspring).
idx_pop <- function(x) {
  if (nchar(x) == 1) grep(x, founders$id)
  else
    match(dat$id[dat$sp == x], as.data.frame(ped)$self)
}

## ----fit1----------------------------------------------------------------
## Avoid estimating BLUPS for which we don't have information
## Otherwise, the run takes much longer (5 hs vs 6 min in this example)

## A[idx_pop('EE'), idx_pop('JJ')]  # This is null: populations are independent
Z_EE <- Z[, idx_pop('EE')]
Z_EJ <- Z[, idx_pop('EJ')]
Z_JJ <- Z[, idx_pop('JJ')]

A_EE <- A[idx_pop('EE'), idx_pop('EE')]
A_EJ <- A[idx_pop('EJ'), idx_pop('EJ')]
A_JJ <- A[idx_pop('JJ'), idx_pop('JJ')]

## Now fit a model with three additive-genetic compnents, 
## by means of generic effects (as only one 'genetic' is allowed in breedR)

res1 <- remlf90(y ~ sp,
                generic = list(
                  E = list(incidence  = Z_EE,
                           covariance = A_EE),
                  J = list(incidence  = Z_JJ,
                           covariance = A_JJ),
                  H = list(incidence  = Z_EJ,
                           covariance = A_EJ)),
                data = dat
)

## ----fit1-summary--------------------------------------------------------
summary(res1)

## ----fit1-predicted-breeding-values--------------------------------------
PBV <- as.matrix(cbind(Z_EE, Z_JJ, Z_EJ)) %*%
  do.call('rbind', lapply(ranef(res1), function(x) cbind(PBV = x, se = attr(x, 'se'))))

ggplot(cbind(dat, PBV), aes(bv, PBV)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = 'darkgray')

## ----fit2----------------------------------------------------------------
## We only want to apply 'dad', 'mum' and 'sca' effects to hybrids,
## and make it zero for non-hybrids. We do so by pre-multiplying by a 
## diagonal indicator matrix
Ind <- diag(dat$sp == 'EJ')

## Build a specific incidence matrices for generic random effects
Z_dad <- Ind %*% as(dat$dad, 'indMatrix')
Z_mum <- Ind %*% as(dat$mum, 'indMatrix')
Z_sca <- Ind %*% as(as.numeric(dat$fam), 'indMatrix')

## The structure variances are diagonal
D <- diag(sum(n.founders))

res2 <- remlf90(y ~ sp,
                generic = list(
                  E = list(incidence  = Z_EE,
                           covariance = A_EE),
                  J = list(incidence  = Z_JJ,
                           covariance = A_JJ),
                  dad = list(incidence = Z_dad,
                             covariance = D),
                  mum = list(incidence = Z_mum,
                             covariance = D),
                  sca = list(incidence = Z_sca,
                             covariance = diag(nlevels(dat$fam)))),
                data = transform(dat)
)

## ----fit2-summary--------------------------------------------------------
summary(res2)

## ----fit2-predicted-breeding-values--------------------------------------
PBV <- as.matrix(cbind(Z_EE, Z_JJ, Z_dad, Z_mum, Z_sca)) %*%
  do.call('rbind', lapply(ranef(res2), function(x) cbind(PBV = x, se = attr(x, 'se'))))

ggplot(cbind(dat, PBV), aes(bv, PBV)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = 'darkgray')

## ----likelihood-profiling------------------------------------------------
## Setup parallel computing
# library(doParallel)
# cl <- makeCluster(2)
# registerDoParallel()
# on.exit(stopCluster(cl))

## Introduce the corresponding scaling factors 
## in the relationship matrix
scale_A <- function(x) {

    ## The pure E subpopulations remains the same
  S <- A
  E.idx <- c(idx_pop('E'), idx_pop('EE'))
  
  ## The pure J subpopulations get multiplied by lambda
  J.idx <- c(idx_pop('J'), idx_pop('JJ'))
  S[J.idx, J.idx] <- A[J.idx, J.idx] * x
  
  ## The hybrids related wuth pure E get a factor of (3+lambda)/4
  S[idx_pop('EJ'), E.idx] <- A[idx_pop('EJ'), E.idx] * (3+x)/4
  S[E.idx, idx_pop('EJ')] <- A[E.idx, idx_pop('EJ')] * (3+x)/4
  
  ## The hybrids related wuth pure J get a factor of (1+3*lambda)/4
  S[idx_pop('EJ'), J.idx] <- A[idx_pop('EJ'), J.idx] * (1+3*x)/4
  S[J.idx, idx_pop('EJ')] <- A[J.idx, idx_pop('EJ')] * (1+3*x)/4

  ## Finally, the hybrids related with other hybrids get a factor of (1+lambda)/2
  S[idx_pop('EJ'), idx_pop('EJ')] <- A[idx_pop('EJ'), idx_pop('EJ')] * (1+x)/2

  return(S)
}

## Condicional likelihood given lambda
cond_lik <- function(x) {
  require(breedR)
  ## Conditional structure matrix
  S <- scale_A(x)
  
  ## Temporarily, let's use only the pure pops
  idx <- c(idx_pop('EE'), idx_pop('JJ'))
  
  suppressWarnings(
    res <- remlf90(y ~ sp,
                   generic = list(
                     E = list(incidence  = Z[dat$sp != 'EJ', idx],
                              covariance = S[idx, idx])),
                   data = dat[dat$sp != 'EJ', ]
    )
  )
  logLik(res)
}

lambda <- seq(.3, 1, length.out = 5)

lik <- sapply(lambda, cond_lik)  # (sequential)
# lik <- foreach(x = seq.int(lambda), .combine = c) %dopar% cond_lik(lambda[x])

ggplot(data.frame(lambda, lik), aes(lambda, lik)) + 
  geom_line()

## ----fit3----------------------------------------------------------------

## Take lambda maximizing the likelihood
lambda0 <- lambda[which.max(lik)]

S <- scale_A(lambda0)

## Temporarily, let's use only the pure pops
idx <- c(idx_pop('EE'), idx_pop('JJ'))

# ## Remove the founders, which I don't want to evaluate
# idx <- -(1:sum(n.founders))

res3 <- remlf90(y ~ sp,
                generic = list(
                  E = list(incidence  = Z[dat$sp != 'EJ', idx],
                           covariance = S[idx, idx])),
                data = dat[dat$sp != 'EJ', ])


## ----fit3-summary--------------------------------------------------------
summary(res3)

## ----fit3-predicted-breeding-values--------------------------------------
PBV <- as.matrix(Z[dat$sp != 'EJ', idx]) %*%
  do.call('rbind', lapply(ranef(res3), function(x) cbind(PBV = x, se = attr(x, 'se'))))

ggplot(cbind(dat[dat$sp != 'EJ', ], PBV), aes(bv, PBV)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = 'darkgray')

