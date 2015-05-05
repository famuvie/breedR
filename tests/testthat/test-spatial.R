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
                           nofill = list(x = 1:6, y = 0:5))),
  # non-consecutive lags
  nclags = list(coord = data.frame(x = c(1,3,5,7,11,13,15), y = c(2,4,8,10,12,14,16)),
              expct = list(fill = list(x = c(1,3,5,7,9,11,13,15), y = c(2,4,6,8,10,12,14,16)),
                           nofill = list(x = c(1,3,5,7,11,13,15), y = c(2,4,8,10,12,14,16)))),
  # non-integer coordinates
  nicoords = list(coord = data.frame(x = c(1.2,2.4,3.6,4.8,7.2), y = c(0.5,1,1.5,2.5,3)),
              expct = list(fill = list(x = c(1.2,2.4,3.6,4.8,6,7.2), y = c(0.5,1,1.5,2,2.5,3)),
                           nofill = list(x = c(1.2,2.4,3.6,4.8,7.2), y = c(0.5,1,1.5,2.5,3)))),
  # negative-valued coordinates 
  nvcoords = list(coord = data.frame( x= c(-9:-5,-3:1), y = c(-11:-9,-7:-1)),
              expct = list(fill = list(x = c(-9:1), y = c(-11:-1)),
                           nofill = list(x = c(-9:-5,-3:1), y = c(-11:-9,-7:-1))))
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

