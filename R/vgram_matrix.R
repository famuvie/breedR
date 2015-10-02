#' Variogram of a matrix
#' 
#' This function is slightly adapted from the homonym function in package \code{fields}.
#' As such, the all the credit is for their authors.
#' 
#' @param dat A matrix.
#' @param R Maximum radius in distance units for computing the variogram.
#' @param dx The spacing between rows of the matrix, in distance units.
#' @param dy The spacing between columns of the matrix, in distance units.
#' 
#' @note 
#' fields, Tools for spatial data
#' Copyright 2004-2013, Institute for Mathematics Applied Geosciences
#' University Corporation for Atmospheric Research
#' Licensed under the GPL -- www.gpl.org/licenses/gpl.html
vgram.matrix <- function(dat, R = 5, dx = 1, dy = 1) {
  # useful function for matching shifted indices
  SI <- function(ntemp, delta) {
    n1 <- 1:ntemp
    n2 <- n1 + delta
    good <- (n2 >= 1) & (n2 <= ntemp)
    cbind(n1[good], n2[good])
  }
  # create all possible separations for a grid up to a distance R
  M <- nrow(dat)
  N <- ncol(dat)
  m <- min(c(round(R/dx), M - 1))
  n <- min(c(round(R/dy), N - 1))
  #
  # all relavent combinations:  note that negative increments are
  # needed as well as  positive ones
  ind <- rbind(as.matrix(expand.grid(0, 1:n)),
               as.matrix(expand.grid(1:m, 0)),
               as.matrix(expand.grid(c(-(m:1), 1:m), 1:n)))
  # distances - only take those within a distance R.
  # and trim everything to this bound
  d <- sqrt((dx * ind[, 1])^2 + (dy * ind[, 2])^2)
  good <- (d > 0) & (d <= R)
  ind <- ind[good, ]
  d <- d[good]
  ind <- ind[order(d), ]
  d <- sort(d)
  #
  # arrays to hold statistics
  nbin <- nrow(ind)
  holdVG <- rep(NA, nbin)
  holdRVG <- rep(NA, nbin)
  holdN <- rep(NA, nbin)
  # loop over each separation
  for (k in 1:nbin) {
    # indices for original and shifted image that are within array bounds
    MM <- SI(M, ind[k, 1])
    NN <- SI(N, ind[k, 2])
    # number of differences and their values
    # holdN[k] <- length(MM) * length(NN)
    BigDiff <- (dat[MM[, 1], NN[, 1]] - dat[MM[, 2], NN[, 
                                                        2]])
    # standard and the  Cressie robust version.
    # modified to handle NAs
    holdN[k] <- sum(!is.na(BigDiff))
    holdVG[k] <- mean(0.5 * (BigDiff)^2, na.rm = TRUE)
    holdRVG[k] <- mean(abs(BigDiff)^0.5, na.rm = TRUE)
  }
  # finish robust estimate Cressie (1993) formula 2.4.12
  holdRVG <- 0.5 * (holdRVG^4)/(0.457 + 0.494 * holdN)
  # collapsed variogram to common distances this what one would look
  # at under the stationary case.
  top <- tapply(holdVG * holdN, d, FUN = "sum", na.rm = TRUE)
  bottom <- tapply(holdN, d, FUN = "sum")
  dcollapsed <- as.numeric(names(bottom))
  vgram <- top/bottom
  #  wipe out pesky row names
  dimnames(vgram) <- NULL
  out <- list(
    vgram = vgram,
    d = dcollapsed,
    n = bottom,
    ind = ind,
    d.full = d,
    vgram.full = holdVG,
    robust.vgram = holdRVG,
    N = holdN,
    dx = dx,
    dy = dy)
  class(out) <- "vgram.matrix"
  return(out)
}
