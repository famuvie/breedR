#' Determine a sensible number of knots
#' 
#' This function computes a reasonable number of inner knots to be used for a 
#' basis of unidimensional B-splines
#' @param n an integer vector of sample sizes
#' @param cutoff a numeric vector of cutoff values
#' @param rate a numeric vector of rates at which the number of knots increases with
#'   the sample size
#' @return An integer vector with the number of knots for each sample size
#' @references Ruppert, D. (2002). Selecting the number of knots for penalized 
#'   splines. \emph{Journal of Computational and Graphical Statistics} 11, 
#'   735–757.
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
#' Given the coordinates of the observations, and the degree, this function puts
#' into place a sensible number of spline knots and computes the incidence 
#' matrix B and the covariance matrix U
#' 
#' @param coord matrix(-like) of observation coordinates
#' @param n.knots numeric. Vector of length two with an integer number of knots 
#'   in each dimension.
#' @param autofill logical. If \code{TRUE} (default) it will try to fill gaps in
#'   the rows or columns. Otherwise, it will treat gaps the same way as adjacent
#'   rows or columns.
#' @param degree integer. Degree of the B-spline polynomials.
build.splines.model <- function (coord, n.knots = NULL, autofill = TRUE, degree = 3) {

  
  # lattice of spatial locations
  # possibly with automatic filling of empty rows or columns
  obs.loc <- loc_grid(coord, autofill)
  
  # Determine the number of (inner) knots for rows and columns
  # or use the number provided by the user
  if(is.null(n.knots))
    n.knots <- determine.n.knots(sapply(obs.loc, length))

  # Spacing between rows/columns
  obs.step <- sapply(obs.loc, function(x) summary(diff(x)))['Median',]
  
  # Place the knots evenly spaced
  # inspired by mgcv::place.knots()
  # but adding a margin before and after the extreme observations
  place.knots <- function(ol, os, nk) {
    if( nk == 1) x <- sum(range(ol))/2
    else x <- seq(head(ol, 1) - os/2, tail(ol, 1) + os/2, length = nk)
    return(x)
  }
  
  # Coordinates of inner knots for rows/columns
  knots.inner <- mapply(place.knots,
                        obs.loc, obs.step, n.knots,
                        SIMPLIFY = FALSE)
  
  # Add n.add additional knots before and after the inner knots
  add.knots <- function(x.inner, n.add, coord) {
    # Use the mean spacing between the inner knots
    # or use the observations range if there are less than two inner knots
    if( length(x.inner) > 1) {
      spacing  <- mean(diff(x.inner))
      extremes <- range(x.inner)
  #     } else if( length(x.inner) == 1) {
  #       spacing  <- diff(range(coord))/2
  #       extremes <- rep(x.inner, 2)
    } else {
      stop('At least two inner knots are needed for each dimension')
    }
    x <- c(extremes[1] - (n.add:1)*spacing,
           x.inner, 
           extremes[2] + (1:n.add)*spacing)
    return(x)
  }

  
  # Return row and column knots in a list
  # as the numbers might be different
  knots <- mapply(add.knots, knots.inner, 3, obs.loc, SIMPLIFY = FALSE)
  
  # Compute incidence matrix B of tensor product of B-spline bases
  # need at least 2*ord -1 knots (typically, 7)
  # but in fact, we need at least 2*ord unless we set outer.ok = TRUE
  # in splineDesign (which we do not want)
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
  
  tensor.sparse <- function (knots, xx, ord) {
    b.x <- as(splineDesign(knots[[1]], xx[, 1], ord = ord), "Matrix")#, sparse=TRUE)
    b.y <- as(splineDesign(knots[[2]], xx[, 2], ord = ord), "Matrix")#, sparse=TRUE)
    # sparse argument was introduced between versions 2.15.0
    # and 3.0.2 of the splines package.
    # sparseness is useful but unnecessary. I prefer to keep this
    # more widely compatible.
    ones.y <- Matrix(1, ncol = ncol(b.y))
    ones.x <- Matrix(1, ncol = ncol(b.x))
    B <- kronecker(b.x, ones.y)*kronecker(ones.x, b.y)
    return(B)
  }
  
  B <- tensor(knots, coord, degree + 1)
  #   B <- tensor.sparse(knots, coord, degree + 1)

  # The sparse incidence matrix weights 3.5 less than the non-sparse version,
  # but it takes longer processing time. Besides, the gmean computation it also
  # takes way longer. And at the end, I need the dense incidence matrix anyway.
  
  U1d <- lapply(sapply(knots, length) - degree - 1, build.splines1d)
  #   U1d <- lapply(sapply(knots, length) - degree - 1, build.splines1d.sparse)
  U <- kronecker(U1d[[1]], U1d[[2]])
  
  #   browser()
  # Scaling so that the characteristic marginal variance equals 1/sigma^2
  # Sorbye and Rue (2014)
  U <- U/gmean(diag(B %*% U %*% t(B)))
  #   U <- U/gmean(Matrix::diag(B %*% U %*% Matrix::t(B)))

  U.sparse <- as(Matrix(tril(U), sparse = TRUE), 'dgTMatrix')
  # Note: The Matrix package counts rows and columns starting from zero
  # Thus, I add 1 to the corresponding columns
  U.values <- cbind(U.sparse@i + 1, U.sparse@j + 1, U.sparse@x)

  # grid and Incidence matrix for plotting purposes
                 
  #   resolution = 51
  #   plot.loc <- sapply(obs.loc, function(x) seq(head(x, 1), tail(x, 1), length=resolution))
  #   plot.grid <- expand.grid(plot.loc[, 1], plot.loc[, 2])
  plot.grid <- expand.grid(obs.loc)
  #   if(!is.null(colnames(coord))) 
  #     colnames(plot.grid) <- colnames(coord)
  plotting <- list(grid = plot.grid,
                   B = tensor.sparse(knots, plot.grid, degree + 1))
  
  return(list(param       = unname(n.knots),
              coord       = coord,
              B           = B,
              U           = U.values,
              Utype       = 'covariance',
              plotting    = plotting))
}


# Build Covariance Matrix in 1D 
build.splines1d <- function(n, model = 'GreenSilverman2003') {
  
  # U matrix (Green & Silverman, 2003)
  temp <- diag(4, n)
  subdiag <- rbind(cbind(0, diag(1, n-1)), 0)
  return((temp + subdiag + t(subdiag))/6)
}


build.splines1d.sparse <- function(n, model = 'GreenSilverman2003') {
  
  # U matrix (Green & Silverman, 2003)
  U <- sparseMatrix(i = c(1:n, 1:(n-1)),
                    j = c(1:n, 2:n),
                    x = c(rep(4, n), rep(1, n-1))/6,
                    dims = c(n, n),
                    symmetric = TRUE)
  
  return(U)
}