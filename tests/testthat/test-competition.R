old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

context("Competition model builder")
########################

data(globulus)
## randomly remove some observations
## and make sure there is some isolated individuals
set.seed(12345)
idx.rm <- sample(nrow(globulus), floor(nrow(globulus)/2))

dat <- globulus[-idx.rm, ]
ped <- globulus[-idx.rm, 1:3]
coord <- globulus[-idx.rm, c('x', 'y')]


gen.spec <- list(model = 'competition',
                 pedigree = ped,
                 id = dat$self,
                 coordinates = coord,
                 competition_decay = 1,
                 pec = TRUE,
                 autofill = TRUE,
                 var.ini = matrix(c(2,-1,-1,2), 2, 2)) 

gen.model <- try(build.genetic.model(gen.spec))
stopifnot(any(apply(gen.model$B[,10:17]>0, 1, sum) == 0))

test_that("The competition model builds without errors", {
  expect_false(inherits(gen.model, "try-error"))
})
