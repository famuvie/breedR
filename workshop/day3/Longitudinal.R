## ----data-preparation, include = FALSE-----------------------------------
data(larix)

if ((current <- packageVersion('breedR')) < "0.10.8") {
  stop(paste('You have breedR v.', current, 'installed.\n',
             'You need the latest version to run this example.\n',
             'Sorry for the inconvenience.'))
}

## Reduce data size for demonstration purposes
dat <- droplevels(larix %>%
                    mutate(yr = ordered(as.numeric(factor(larix$map)))) %>% 
                    filter(rep %in% c("1", "2", "7"),
                           yr %in% c("4", "5", "7", "8", "10", "11", "13", "14", "15")))


## ----data-description, echo = FALSE--------------------------------------
with(dat, table(yr, rep))


## ----response-by-mean-annual-precipitation, echo = FALSE, fig.width=8, fig.height = 5----
ggplot(transform(dat, map = factor(map)), aes(x, y)) +
  geom_tile(aes(fill = LAS)) +
  facet_wrap(~ map)

## ----mean-response-by-year, echo = FALSE---------------------------------
larix %>%
  group_by(yr, rep) %>%
  summarise(MLAS = mean(LAS)) %>%
  ggplot(aes(yr, MLAS, group = rep)) +
  geom_line()


## ----computation-legendre, echo = FALSE----------------------------------

#' Evaluate Legendre polynomials
#' 
#' Evaluate the Legendre polynomials up to certain order on a given set of 
#' values. Normalization is done within the function.
#' 
#' By default,
#' 
#' @param x numeric. Vector of values for evaluation after normalization into 
#'   [-1, 1]
#' @param ord integer. Maximum order of the polynomials to evaluate.
#' @param scale vector of length 2. Scale to use for normalization. See Details.
#'   
#'   By default, the empirical range of \code{x} is linearly scaled into 
#'   [-1, 1]. However, sometimes it is useful to give a predefined range. For
#'   example, when the legendre values are required at several different sets of
#'   points x but in the same consistent scale. Note that \code{range = 0:1}
#'   implies no re-scaling at all.
#' @export
legendre_values <- function (x, ord, scale = range(x)) {
  # Normalize x
  if( diff(scale) > 0 ) {
    w <- 2*(x - scale[1])/diff(scale) - 1
  } else {
    # constant vector or one value only
    w <- 0
  }
  # Legendre polynomials list up to order ord
  p.list <- legendre.polynomials(ord, normalized = TRUE)
  
  # Polynomial values for each order
  v.list <- polynomial.values(p.list, w)
  
  # Return matrix-wise
  return(do.call(cbind, v.list))
}


## ----legendre-polynomials, echo = FALSE----------------------------------
# time points 
x <- seq(from = -1, to = 1, length = 101)
# basis matrix evaluated at the time points t
# an intercept column is included
Phi <- legendre_values(x, ord = 3)

# # check the orthogonality
# crossprod(Phi) # is equal to I_4 up to rounding errors

# plot the basis system excluding the constant column
datPhi <- data.frame(x = x, Phi)
colnames(datPhi) <- c('x', paste('order', 0:3, sep = '_'))
datPhi <-
  datPhi %>% 
  tidyr::gather('order', 'y', 2:5) %>% 
  mutate(order = as.factor(extract_numeric(order)))


ggplot(datPhi, aes(x, y, color = order)) +
  geom_line()

## ----breedR-hack, echo = FALSE, warning = FALSE, message = FALSE---------
## Order of polynomials and longitudinal variable
ord <- 4
long_var <- dat$map

# Start by fitting the model as if we had one observation per individual
# (i.e. as if it wasn't longitudinal)
# This is only in order to get a first sketch of the model structure

res_dummy <- remlf90(LAS ~ 1 + rep + yr,
                     random = ~ bl,
                     genetic = list(model = 'add_animal',
                                    pedigree = dat[, 1:3],
                                    id = 'self'),
                     data = dat)


# Now we are going to "manually" modify the relevant internal structures
eff_dummy <- res_dummy$effects

eff <- within(eff_dummy, {
  # The 'genetic' effect will no longer be one value, but (ord + 1)
  genetic$effects <- rep(genetic$effects, ord + 1)
  names(genetic$effects) <- paste('a', seq_len(ord + 1), sep = '_')
  genetic$cov.ini <- diag(as.vector(genetic$cov.ini), ord + 1)
#   genetic$cov.ini <- matrix(as.vector(genetic$cov.ini)/2, ord+1, ord+1)
#   diag(genetic$cov.ini) <- 2*genetic$cov.ini[1,1]
})

pf90_dummy <- breedR:::progsf90(res_dummy$mf, eff, opt = c("sol se"), res.var.ini = 1)
class(pf90_dummy) <- c(class(pf90_dummy), 'list')
rendered <- breedR:::renderpf90.breedr_modelframe(eff, ntraits = 1)

pf90 <- within(pf90_dummy, {
  # ord + 1 covariates nested in animal
  # new position for the animal id
  parameter$effects$genetic <- 
    paste0(gsub('cross', 'cov', parameter$effects$genetic),
         max(rendered$genetic$pos) + 1)
  
  # include legendre coefficients in data
  data[, rendered$genetic$pos] <-
    legendre_values(long_var, ord)
  data <- cbind(data, rendered$genetic$data[,1])
})

# Write files and run EMREML
tmpdir <- tempdir()
breedR:::write.progsf90(pf90, dir = tmpdir)

cdir <- setwd(tmpdir)

reml.out <- system2(file.path(breedR.getOption("breedR.bin"), 'remlf90'), 
                    input  = 'parameters',
                    stdout = TRUE)

# Parse solutions
res <- breedR:::parse_results('solutions', eff, res_dummy$mf, reml.out, 'em', res_dummy$call)
class(res) <- c('breedR', 'remlf90')

setwd(cdir)


## ---- echo = FALSE, warning = FALSE--------------------------------------

# predicted coefficients of breeding values
# in a continuous scale
lims <- range(long_var)
x <- seq(lims[1], lims[2], length = 11)
Zmat <- legendre_values(seq(lims[1], lims[2], length = 11), ord)
coefs <- sapply(res$ranef[-1], function(x) x$value)
pred_BV <- coefs %*% t(Zmat)
colnames(pred_BV) <- x
pred_BV <- data.frame(individual = seq_len(nrow(pred_BV)), pred_BV)
pred_BV <- pred_BV %>% 
  gather('precip', 'value', -1) %>% 
  mutate(precip = extract_numeric(precip))


inc_mat <- legendre_values(dat$mat, 4)
ped <- build_pedigree(1:3, data = dat)
Amat <- pedigreemm::getA(ped)
map <- attr(ped, 'map')
ind.codes <- map[pred_BV$individual]

# dat %>%
#   select(self:mum) %>% 
#   unique() %>% 
#   mutate(fam = factor(dad:mum))


## Find most interactive individuals
# plot(t(coefs))
# apply(t(coefs), 2, summary)

## Identify individuals that maximise each of the coefficients
fam.idx <- apply(coefs, 2, which.max)

ggplot(pred_BV, aes(precip, value, group = individual)) +
  geom_line() + 
  geom_line(data = pred_BV[pred_BV$individual %in% fam.idx, ], color = 'red', lwd = 1)



# true vs. predicted coefficients
# ggplot(data.frame(gen = factor(rep(c('founder', 'progeny'),
#                                    times = c(sum(N$fnd), N$ind))),
#                   coef = factor(rep(1:(ord+1), each = sum(N$fnd) +N$ind)),
#                   true_BV = as.vector(true_BV),
#                   pred_BV = as.vector(pred_BV)),
#        aes(true_BV, pred_BV)) + 
#   geom_point(aes(col = coef)) + 
#   geom_abline(int = 0, sl = 1, col = 'darkgray') +
#   coord_fixed()
# 
# 
# 
# # fitted vs. observations
# ggplot(data.frame(fit = fitted(res),
#                   obs = dat$phe),
#        aes(fit, obs)) + 
#   geom_point() + 
#   geom_abline(int = 0, sl = 1, col = 'darkgray') + 
#   coord_fixed()
# 
# 
# # variance components
# ggplot(data.frame(name = c(paste('a', 1:(ord+1), sep ='_'), 'resid'),
#                   est = with(res$var, c(diag(genetic), res=Residual)),
#                   true = c(rep(sigma2$a, ord+1), sigma2$e)),
#        aes(true, est)) + 
#   geom_point(aes(col = name)) + 
#   geom_abline(int = 0, sl = 1, col = 'darkgray') +
#   coord_fixed()
# 
# 
# # Some BVs in time
# # The points represent the times where the phenotype has been measured
# # for each individual
# n_show = 3
# idx <- sample(sum(N$fnd)+1:N$ind, n_show)
# x <- seq(t_range[1], t_range[2], length = 101)
# true_BVf <- legendre_values(x, ord) %*% t(true_BV[idx,])
# pred_BVf <- legendre_values(x, ord) %*% t(pred_BV[idx,])
# BVf <- transform(x,
#                  t = x,
#                  ind = factor(rep(idx, each = length(x), times = 2)),
#                  type = factor(rep(c('true', 'pred'), each = length(x)*n_show)),
#                  BV = c(as.vector(true_BVf), as.vector(pred_BVf)))
# 
# BVt <- mapply(function(i, t) legendre_values(t, ord,
#                                              scale = t_range) %*% pred_BV[i,],
#               idx,
#               tobs[idx],
#               SIMPLIFY = FALSE)
# 
# BV_pt <- do.call(rbind,
#                  lapply(1:n_show,
#                         function(i) data.frame(ind = factor(idx[i]),
#                                                t = tobs[[idx[i]]],
#                                                BV = BVt[[i]])))
# 
# ggplot(BVf, aes(t, BV)) + 
#   geom_line(aes(col = ind, lty = type)) + 
#   geom_point(aes(col = ind), data = BV_pt)



