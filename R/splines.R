#' Determine a sensible number of knots
#' 
#' This function computes a reasonable number of inner knots
#' to be used for a basis of unidimensional B-splines
#' @param n An integer vector of sample sizes
#' @return An integer vector with the number of knots for each sample size
#' @references Ruppert, D. (2002). Selecting the number of knots for penalized splines. \emph{Journal of Computational and Graphical Statistics} 11, 735–757.
determine.n.knots <- function(n, cutoff = 4, rate = 0.3) {
  # # Inspired by hisemi::n.knots
  # as.integer(trunc(pmin(n, cutoff + pmax(0, n - cutoff)^rate)))
  # but taking into account that we will add three more knots at each side
  if(any(n < 7)) stop('Too few data for a bidimensional spatial model')
  nk <- trunc(pmin(n - 6, cutoff + pmax(0, n - 6 - cutoff)^rate))
  return(nk)
}


#' Build splines model
#' 
#' Given the coordinates of the observations, and the degree,
#' this function puts into place a sensible number of spline knots
#' and computes the incidence matrix B and the covariance matrix U
build.splines.model <- function (coord, degree = 3) {
  # Determine the number of (inner) knots for rows and columns
  # TODO: Let the user fix these numbers
  obs.loc <- lapply(as.data.frame(coord), function(x) sort(unique(x)))
  n.knots <- determine.n.knots(sapply(obs.loc, length))
  obs.step <- sapply(obs.loc, function(x) summary(diff(x)))['Median',]
  
  # Place the knots evenly spaced
  # inspired by mgcv::place.knots()
  # but adding a margin before and after the extreme observations
  knot.step <- mapply(function(ol, os, nk) diff(range(ol) + c(-1, 1)*os/2)/(nk-1), obs.loc, obs.step, n.knots)
  knots.inner <- mapply(function(ol, os, nk) seq(head(ol, 1) - os/2, tail(ol, 1) + os/2, length = nk), obs.loc, obs.step, n.knots, SIMPLIFY = FALSE)
  
  # Add three additional knots before and after
  add.knots <- function(x.inner, n.add) {
    # Use the mean spacing between the inner knots
    spacing <- mean(diff(x.inner))
    x <- c(x.inner[1] - (n.add:1)*spacing, x.inner, tail(x.inner, 1) + (1:n.add)*spacing)
    return(x)
  }
  knots <- lapply(knots.inner, add.knots, 3)
  
  # Compute incidence matrix B of tensor product of B-spline bases
  tensor <- function (knots, xx, ord) {
    b.x <- splineDesign(knots[[1]], xx[, 1], ord = ord)#, sparse=TRUE)
    b.y <- splineDesign(knots[[2]], xx[, 2], ord = ord)#, sparse=TRUE)
      # sparse argument was introduced between versions 2.15.0
      # and 3.0.2 of the splines package.
      # sparseness is useful but unnecessary. I prefer to keep this
      # more widely compatible.
    ones.y <- matrix(1, ncol = ncol(b.y))
    ones.x <- matrix(1, ncol = ncol(b.x))
    B <- kronecker(b.x, ones.y)*kronecker(ones.x, b.y)
    return(B)
  }
  
  B <- tensor(knots, coord, degree + 1)
  
  # Compute U matrix (Green & Silverman, 2003)
  Sigma.marginal <- function(n) {
    temp <- diag(4, n)
    subdiag <- rbind(cbind(0, diag(1, n-1)), 0)
    return((temp + subdiag + t(subdiag))/6)
  }
  U <- kronecker(Sigma.marginal(length(knots[[1]])-degree-1), 
                 Sigma.marginal(length(knots[[2]])-degree-1))
  
  U.sparse <- as(Matrix(tril(U), sparse = TRUE), 'dgTMatrix')
  # Note: The Matrix package counts rows and columns starting from zero
  # Thus, I add 1 to the corresponding columns
  U.values <- cbind(U.sparse@i + 1, U.sparse@j + 1, U.sparse@x)

  
  # grid and Incidence matrix for plotting purposes
                 
  resolution = 101
  plot.loc <- sapply(obs.loc, function(x) seq(head(x, 1), tail(x, 1), length=resolution))
  plot.grid <- expand.grid(plot.loc[, 1], plot.loc[, 2])
  if(!is.null(colnames(coord))) 
    colnames(plot.grid) <- colnames(coord)
  plotting <- list(grid = plot.grid,
                   B = tensor(knots, plot.grid, degree + 1))
  
  return(list(B = B, U = U.values, plotting = plotting))  
}
