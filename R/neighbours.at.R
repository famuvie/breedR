#' @describeIn neighbours.at Returns the matrix of neighbours in the specified
#'   direction
#' @export
neighbours.at.matrix <-function(x, dir) {
  stopifnot(is.character(dir))
  if( length(dir) == 1 ){
    dimx <- dim(x)
    switch(EXPR = dir,
           N = rbind(NA, x[-dimx[1], ]),
           S = rbind(x[-1, ], NA),
           E = cbind(x[, -1], NA),
           W = cbind(NA, x[, -dimx[2]]),
           NE = cbind(rbind(NA, x[-dimx[1], -1]), NA),
           SE = cbind(rbind(x[-1, -1], NA), NA),
           SW = cbind(NA, rbind(x[-1, -dimx[2]], NA)),
           NW = cbind(NA, rbind(NA, x[-dimx[1], -dimx[2]]))
    )
  } else {
    sapply(dir, function(d) neighbours.at(x, d))
  }
}

#' @describeIn neighbours.at Returns the list of matrices of neighbours in the
#'   specified direction
#' @export
neighbours.at.list <- function(x, dir) {
  # Check that it is a list of matrices
  stopifnot(all(lapply(x, class) == 'matrix'))
  lapply(x, neighbours.at.matrix, dir)
}
