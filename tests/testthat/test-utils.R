### Test the auxiliar functions in utils.R ###

old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

context("Auxiliar functions")

## building matrices by blocks

test_that('Diagonal-binding function dbind() works as expected', {
  x <- list(matrix(0+1:6, 2,3), matrix(10+1:9, 3,3), 20+1:4)
  fillin <- 0
  res <- dbind(x, fillin)
  
  expect_identical(res[1:2, 1:3], x[[1]])
  expect_identical(res[3:5, 4:6], x[[2]])
  expect_identical(res[6:9, 7], x[[3]])
  expect_identical(res[3, 1], fillin)
})
