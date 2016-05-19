
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
    res <- try(renderpf90.matrix(x$m))
    
    expect_true(!(failed <- inherits(res, 'try-error')))
    
    if (!failed) {
      expect_equal(res, x$r)
    }
  }

})


test_that('renderpf90.breedr_modelframe() renders a single-trait breedr_modelframe correctly', {

  # TODO...  
  # testdat <- transform(
  #   expand.grid(x = 1:4, y = 1:4, KEEP.OUT.ATTRS = FALSE),
  #   z = rnorm(16),
  #   mu = 1)
  # 
  # bc <- call('remlf90', fixed = phe_X~1, data = quote(as.data.frame(m1)))
  # str(build.mf(bc))
  # 
  # 
  # breedrmf <- build.effects(mf = build.mf(bc),
  #                           genetic = NULL,
  #                           spatial = NULL,
  #                           generic = NULL)
  # renderpf90.breedr_modelframe(breedrmf, ntraits = 1)
})
