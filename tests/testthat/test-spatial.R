old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

### Internal spatial functions ###
context("Internal spatial functions")

## list of test coordinates with expected result
tcl <- list(
  # holes in rows and columns
  holes = list(coord = data.frame(x = c(1:3, 5:8), y = c(1:2, 5:9)),
               expct = list(fill = list(x = 1:8, y = 1:9),
                            nofill = list(x = c(1:3, 5:8), y = c(1:2, 5:9)))),
  # replicates in same position
  reps = list(coord = data.frame(x = c(1:3, 3:6), y = c(0:4, 4:5)),
              expct = list(fill = list(x = 1:6, y = 0:5),
                           nofill = list(x = 1:6, y = 0:5)))
  
  # TODO: 
  #   - non-consecutive lags
  #   - non-integer coordinates
  #   - negative-valued coordinates
  )

## Debug: plot datasets
# for(x in tcl) 
#   print(ggplot(x$coord, aes(x, y)) + geom_point())

  
test_that("loc_grid behaves as expected", {
  for(x in tcl) {
    expect_equal(loc_grid(x$coord, autofill = TRUE),
                 x$expct$fill)
    expect_equal(loc_grid(x$coord, autofill = FALSE),
                 x$expct$nofill)
  }
})

## TODO:
##   - test functions fill_holes() and fill_hole() (see defs in ../R/spatial.R)