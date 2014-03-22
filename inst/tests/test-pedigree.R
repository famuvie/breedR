#### pedigree building and checking ####
context("Pedigree")

# Use the pedigree in data(m4) and shuffle the codes
data(m4)
ped <- as.data.frame(m4)[, c('self', 'dad', 'mum')]
test_that('The pedigree from m4 is not complete, but otherwise correct', {
  expect_that(!check_pedigree(ped)['full_ped'], is_true());
  expect_that(all(check_pedigree(ped)[-1]), is_true())
})

# Generate a crazy map
mcode <- max(ped, na.rm = TRUE)
map <- rep(NA, mcode)
set.seed(1234)
map <- sample(10*mcode, size = mcode)

# Generate a crazy pedigree that fails all checks
ped_shuffled <- sapply(ped, function(x) map[x])
# Introduce some unknown parents either with NA or with 0
ped_shuffled[, 2:3][sample(2*nrow(ped), 200)] <- c(0, NA)

test_that('The shuffled pedigree fails all checks', {
  expect_that(all(!check_pedigree(ped_shuffled)), is_true())
})

# Reorder and recode 
ped_fix <- suppressWarnings(build_pedigree(1:3, data = ped_shuffled))
test_that('build_pedigree() fixes everything', {
  expect_that(all(check_pedigree(ped_fix)), is_true())
})


# Check that remlf90 handles correctly recoded pedigrees
# by comparing the genetics evaluations of a dataset with or without
# a shuffled pedigree

data(m1)
dat <- as.data.frame(m1)
ped <- get_pedigree(m1)

res_ok <- try(
  remlf90(fixed = phe_X ~ sex, 
          genetic = list(model = 'add_animal', 
                         pedigree = ped,
                         id = 'self'), 
          data = dat),
  silent = TRUE
)

# Shuffle the pedigree
mcode <- max(as.data.frame(ped), na.rm = TRUE)
map <- rep(NA, mcode)
set.seed(1234)
map <- sample(10*mcode, size = mcode)
m1_shuffled <- m1
m1_shuffled$Data[, 1:3] <- sapply(as.data.frame(ped), function(x) map[x])

ped_fix <- suppressWarnings(
  build_pedigree(1:3, data = as.data.frame(get_pedigree(m1_shuffled)))
)

res_shuffled <- try(
  remlf90(fixed = phe_X ~ sex,
          genetic = list(model = 'add_animal', 
                         pedigree = ped_fix,
                         id = 'self'), 
          data = as.data.frame(m1_shuffled)),
  silent = TRUE
)

# Except the call, and the reml output everything must be the same
# Update: also need to omit the shuffled random effects estimations
# which should be the same, but reordered
test_that('remlf90 handles recoded pedigrees correctly', {
  omit.idx <- match(c('call', 'reml', 'ranef'), names(res_ok))
  expect_that(res_ok[-omit.idx], equals(res_shuffled[-omit.idx]))
})