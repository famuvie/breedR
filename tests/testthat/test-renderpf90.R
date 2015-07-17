old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

### Test the rendering to pf90 format ###

context("Render PF90")

test_that('renderpf90.matrix() renders different matrix types', {

  mat.list <- list(
    ## one effective column, some rows full-zero
    list(m = diag(c(2, 0, 4)),
         r = matrix(c(2, 0, 4, 1, 0, 3), ncol = 2))
    ,
    ## two effective columns, some rows full-zero
    ## some column fill-in with zero also needed
    list(m = rbind(c(11, 13, 0), 0, c(0, 0, 14)),
         r = rbind(c(11, 13, 1, 2), 0, c(14, 0, 3, 0)))
  )
  
  for (x in mat.list) {
    res <- try(renderpf90.matrix(x$m), silent = TRUE)
    
    expect_true(!(failed <- inherits(res, 'try-error')))
    
    if (!failed) {
      expect_equal(res, x$r)
    }
  }

})

