## test data directory
testdata <- system.file("testdata", package = "breedR")

load_res <- function(key, dir = testdata) {
  fn <- paste0("res_", key, ".rds")
  readRDS(file.path(dir, fn))
}