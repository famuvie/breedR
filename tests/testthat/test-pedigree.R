#### pedigree building and checking ####

context("Pedigree infrastructure")

# Retrieve pedigree from remlf90 objects

test_that('get_pedigree() returns NULL when there is no genetic effect', {
  res <- load_res("fixonly")
  expect_true(is.null(get_pedigree(res)))
})


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

