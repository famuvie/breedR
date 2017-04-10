# Generate (or update) testdata used by unit tests
message("Generating test data ...")

## target directory
testdata <- system.file(file.path("inst", "testdata"), package = "breedR")


#### fitted models ####

res  <- suppressMessages(
  list(
    ## A simple model with one fixed effect
    fixonly = remlf90(
      fixed  = phe_X ~ gg,
      data = globulus
    ),
    ## A spatial blocks model
    blk = remlf90(
      fixed  = phe_X ~ 1,
      spatial = list(
        model = 'blocks', 
        coord = globulus[, c('x','y')],
        id = 'bl'
      ),
      data = globulus
    ),
    ## An spatial autoregressive model with one random effect
    ar = remlf90(
      fixed  = phe_X ~ 1,
      random = ~ gg,
      spatial = list(
        model = 'AR', 
        coord = globulus[, c('x','y')],
        rho = c(.85, .8)
      ), 
      data = globulus
    ),
    ## An spatial splines (2x2 knots) model fitted with EM
    spl = remlf90(
      fixed  = phe_X ~ 1,
      spatial = list(
        model = 'splines', 
        coord = globulus[, c('x','y')], 
        n.knots = c(2, 2)
      ), 
      data = globulus,
      method = 'em'
    ),
    ## A genetic-AR model with a fixed effect
    ped_ar = remlf90(
      fixed  = phe_X ~ gg,
      genetic = list(
        model = 'add_animal', 
        pedigree = globulus[,1:3],
        id = 'self'
      ), 
      spatial = list(
        model = 'AR', 
        coord = globulus[, c('x','y')],
        rho = c(.85, .8)
      ), 
      data = globulus
    )
  )
)


for (idx in seq_along(res)){
  fn <- paste0("res_", names(res)[idx], ".rds")
  saveRDS(res[[idx]], file = file.path(testdata, fn))
}

