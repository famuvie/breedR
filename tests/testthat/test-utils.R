### Test the utility functions ###

context("matrix.short16")


## Test matrix.short16
M1<-matrix(c(rep(3,4),rep(1,4),1,2,4,6,8,9,11,12),nrow=1,ncol=16,byrow=TRUE)
M2<-matrix(c(1,5,-6,4,5,3,-1,8,1,3,7,5,9,11,15,6,0,0,rep(1,5),0,1,2,4,6,8,9,11,12,0,rep(2,6),0,2,3,4,6,7,8,11,12),
           nrow=3,ncol=16,byrow=TRUE)
Mtest1<-matrix.short16(M1)
Mtest2<-matrix.short16(M2)

test_that("The function matrix.short16 behaves as expected", {
  i=1
  j=3
  expect_equal(M1[i,j],Mtest1[i,M1[i,j+8]])
  i=2
  j=5
  expect_equal(M2[i,j],Mtest2[i,M2[i,j+8]])
  i=3
  j=5
  expect_equal(M2[i,j],Mtest2[i,M2[i,j+8]])
})

