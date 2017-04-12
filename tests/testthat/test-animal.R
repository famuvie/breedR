
data(m1)
dat <- as.data.frame(m1)
ped <- get_pedigree(m1)

#### Context: Building additive-genetic models ####
context("Building additive-genetic models")

test_that("Correctly builds structure of additive_genetic_animal component", {
  dat <- data.frame(id = 1:4,
                    sire = c(11, 11, 2, 3),
                    dam  = c(12, NA, 1, 12))
  ## Recoded pedigree with further unobserved descendants
  ## (with will come later in the codification)
  ped <- suppressWarnings(build_pedigree(1:3, data = rbind(dat, c(5, 1, 2))))
  aga <- try(breedR:::additive_genetic_animal(ped, dat$id))
  
  expect_true(!inherits(aga, "try-error"))
  expect_is(aga, 
            c("additive_genetic_animal", 
              "additive_genetic",
              "genetic",
              "random", 
              "breedr_effect"))
  expect_named(aga, 
               c("incidence.matrix", 
                 "structure.matrix",
                 "structure.type", 
                 "pedigree"))
  expect_identical(dim(aga$incidence.matrix),
                   c(nrow(dat), nrow(as.data.frame(ped))))
  expect_identical(dim(aga$structure.matrix),
                   rep(nrow(as.data.frame(ped)), 2))
  expect_identical(aga$structure.type, 'covariance')
  expect_identical(aga$pedigree, ped)
})
