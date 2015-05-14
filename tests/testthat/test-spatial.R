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

## Tests for the build_grid function
x1 = c(-1:3, 5:8)
y1 = c(-5:-2, 5:9)
A<-matrix(c(x1,y1),ncol=2)
x2 = c(2.1,4.2,6.3,10.5,14.7,8.4,16.8)
y2 = c(1,4,13,16,19,22,7)
B<-matrix(c(x2,y2),ncol=2)
C<-matrix(c(6:18),ncol=1)
D<-matrix(c(1:24),ncol=4)
E<-matrix(c(1,3,4,8,11,13,15,1,4,3,5,11,7,6),ncol=2)

## TODO:
## There are a couple of failing tests in what follows (check)
## Particularly, loc_grid(B, autofill = TRUE) returns some duplicated
## coordinates in the first dimension. 
## Check why is this happening and fix.

test_that("The function stops if the number of columns does not equal 2",{  
expect_error(build_grid(C))
expect_error(build_grid(D))
})
  
test_that("The function stops if the grid is not regular",{  
  expect_error(build_grid(E))
})
  
test_that("The function returns the right origin",{
  expect_true(all(build_grid(A)$origin==c(min(A[,1]),min(A[,2]))))
  expect_true(all(build_grid(B)$origin==c(min(B[,1]),min(B[,2]))))
})
  
test_that("The function returns the right step",{
    expect_true(all(build_grid(A)$step==c(1,1)))
    expect_true(all(build_grid(B)$step==c(2.1,3)))
})
  
test_that("The function returns the right length",{
    expect_true(all(build_grid(A)$length==c(10,15)))
    expect_true(all(build_grid(B)$length==c(8,8)))
})
  
test_that("The function returns the right idx",{
    expect_true(all(build_grid(A)$idx==c(1,12,23,34,105,117,128,139,150)))
    expect_true(all(build_grid(B)$idx==c(1,10,35,45,55,60,24)))
})
  
  