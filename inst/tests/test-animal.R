data(m1)
dat <- as.data.frame(m1)
ped <- get_pedigree(m1)

 res <-  remlf90(fixed = phe_X ~ sex, 
          genetic = list(model = 'add_animal', 
                         var.ini = 10, 
                         pedigree = ped,
                         id = 'self'), 
          data = dat,
          method = 'ai')

#### Context: Animal Models ####
context("Animal Models") 

fixed_models <- list(phe_X ~ sex)

# Run REML and lm and save estimates and MLEs
run_model <- function(m, data = dat, method) {
  res.reml <- try(
    remlf90(fixed = m,
            genetic = list(model = 'add_animal', 
                           var.ini = 10, 
                           pedigree = ped,
                           id = 'self'), 
            data = data,
            method = method),
    silent = TRUE)
  return(res.reml)
}

# Compare progsf90 and lm results
run_expectations <- function(m, data = dat, method) {
  res <- run_model(m, data, method)
  
  # It runs without errors 
  test_that("The animal model runs without errors", {
    expect_that(!inherits(res, "try-error"), is_true())
  })
  
  # TODO:
  # other checks, like:
  #  - compare the estimated genetic and residual vaiances with true values
  #  - compare the estimated and true Breeding Values
  
}



# Run expectations for all models and methods
test_that("remlf90() estimates matches lm()'s", {
  lapply(fixed_models, run_expectations, method = 'em')
})

test_that("airemlf90() estimates matches lm()'s", {
  lapply(fixed_models, run_expectations, method = 'ai')
})

