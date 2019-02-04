## test data directory
testdata <- system.file("testdata", package = "breedR")

load_res <- function(key, dir = testdata) {
  fn <- paste0("res_", key, ".rds")
  readRDS(file.path(dir, fn))
}


## Generate sample REML output files

## 10-variate model.
## large residual covariance matrix (10x10) which wraps lines
if (!file.exists(tf <- file.path(testdata, "airemlf90_log_4.txt"))) {
  
  tdat <- data.frame(replicate(n = 10, rnorm(1e2)))
  res <- remlf90(
    cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) ~ 1,
    data = tdat,
    method = "ai"
  )
  writeLines(res$reml$output, tf)
}
