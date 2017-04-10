#### pedigree building and checking ####

context("Pedigree")

# Toy dataset with silly pedigree
test.dat <- data.frame(matrix(sample(100, 15), 5, 3,
                              dimnames = list(NULL, c('self', 'sire', 'dam'))),
                       y = rnorm(5))
ped.fix <- suppressWarnings(build_pedigree(1:3, data = test.dat))
test.res <- try(
  suppressMessages(
    suppressWarnings(
      remlf90(y~1,
              genetic = list(model = 'add_animal',
                             pedigree = test.dat[, 1:3],
                             id = 'self'),
              data = test.dat)
    )
  ),
  silent = TRUE
)

test_that('remlf90() builds and recodes the pedigree', {
  expect_false(inherits(test.res, 'try-error'))
})

test_that('get_pedigree() returns the recoded pedigree', {
  expect_identical(ped.fix, get_pedigree(test.res))
})


# Check that remlf90 handles correctly recoded pedigrees
# by comparing the genetics evaluations of a dataset with or without
# a shuffled pedigree

data(m1)
dat <- as.data.frame(m1)
ped <- get_pedigree(m1)

res_ok <- try(
  suppressMessages(
    remlf90(fixed = phe_X ~ sex, 
            genetic = list(model = 'add_animal', 
                           pedigree = ped,
                           id = 'self'), 
            data = dat)
  )
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
  suppressMessages(
    remlf90(fixed = phe_X ~ sex,
            genetic = list(model = 'add_animal', 
                           pedigree = ped_fix,
                           id = 'self'), 
            data = as.data.frame(m1_shuffled))
  )
)

# Except the call, and the reml output everything must be the same
# Update: also need to omit the shuffled random effects estimations
# which should be the same, but reordered
test_that('remlf90 handles recoded pedigrees correctly', {
  omit.idx <- match(c('call', 'effects', 'reml', 'ranef'), names(res_ok))
  expect_that(res_ok[-omit.idx], equals(res_shuffled[-omit.idx]))
})
