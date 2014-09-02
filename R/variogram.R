# TODO
# - number of observations per variogram cell 
# - restrict R to cells where N > 30

#' Empirical variograms of residuals
#' 
#' Computes isotropic and anisotropic empirical variograms from the residuals of
#' a breedR model.
#' 
#' An empirical variogram computes the mean squared differences between 
#' observations separated by a vector \code{h}, for all possible \code{h}. An 
#' \code{isotropic} variogram assumes that the underlying process depends only 
#' on the (euclidean) \emph{distance} between points, but not on the orientation
#' or direction of \code{h}. At the other end, an \code{anisotropic} variogram 
#' assumes that the process might depend on the orientation, but not on the 
#' direction. Finally, \code{heat} and \code{perspective} are different 
#' representations of a variogram which assumes that the process depends only on
#' the absolute distance between rows and columns.
#' 
#' Unless \code{coord} or \code{z} are specified by the user, \code{variogram} 
#' builds the variogram with the residuals of the model fit in \code{x}. If
#' \code{coord} or \code{z} are specified, then the spatial coordinates or the
#' residuals are respectively overrided.
#' 
#' This function assumes that there is at most one observation per spatial 
#' location. Otherwise, are observations measured at different times? should a 
#' spatial-temporal variogram be fitted?
#' 
#' @param x a \code{breedR} result.
#' @param plot character. What type of variogram is to be plotted. The default 
#'   is 'all'. Other options are 'isotropic', 'heat', 'perspective' and 
#'   'anisotropic'. See Details.
#' @param R numeric. Radius of the variogram
#' @param coord (optional) a two-column matrix with coordinates of observations
#' @param z (optional) a numeric vector of values to be represented spatially
#'   
#'   
#' @importFrom fields vgram.matrix
#' @export
variogram <- function(x, plot = c('all', 'isotropic', 'anisotropic', 'perspective', 'heat', 'none'), R, coord, z) {
  
  if( missing(x) ) {
    if( (missing(coord) | missing(z)) ) {
      stop(paste('This function works on models fitted with breedR,',
                 'or with both arguments coord and z.\n',
                 sep='\n'))
    } else {
      coord <- as.matrix(coord)
      stopifnot(is.numeric(coord))
      stopifnot(ncol(coord) == 2)
      stopifnot(nrow(coord) > 1)
    }
  } else {
    stopifnot(inherits(x, 'breedR'))
  }
  
  plot <- match.arg(plot)
  # We need to place the residuals in matrix form
  # where rows and columns have the same spacing
  
  # The user can provide coordinates for non-spatial remlf90 models
  # Otherwise, use coordinates of fitted dataset
  if(missing(coord)) {
    # Coordinates in distance units (e.g. m)
    coord <- coordinates(x)
  }
  
  # The user can provide a custom vector to be plotted
  # Otherwise, plot variogram of residuals
  if(missing(z)) {
    z <- residuals(x)
  }
  
  # Check: Only one observation per spatial unit
  if( any(duplicated(coord)) )
    stop('More than one observation detected in at least one spatial unit.\n')
  
  # A regular grid containind all the location units
  loc.grid <- loc_grid(coord, autofill = TRUE)
  
  # Spacing between rows/columns (in distance units)
  step <- sapply(loc.grid, function(x) summary(diff(x)))['Median',]
  
  # Coordinates in row/col numbers
  coord1 <- as.data.frame(mapply(factor,
                                 x = as.data.frame(coord),
                                 levels = loc.grid,
                                 SIMPLIFY = FALSE))
  
  # Matrix of residuals by rows and cols
  dat <- Matrix::sparseMatrix(i = as.integer(coord1[[1]]),
                              j = as.integer(coord1[[2]]),
                              x = z
                              , dimnames = lapply(coord1, levels)
  )
  
  # Maximum radius for the variogram (in distance units)
  if( missing(R) ) R <- max(coord/sqrt(2))
  
  # Variogram computation
  out <- vgram.matrix(as.matrix(dat), R = R, dx = step[1], dy = step[2])

  # Anisotropic variogram
  #   plot(out)  # from package fields
  # TODO: change gradient colours.
  dat.aniso <- with(out, data.frame(x = ind[, 1] * dx,
                                    y = ind[, 2] * dy,
                                    z = vgram.full))

  # Isotropic variogram
  dat.iso <- data.frame(distance = out$d,
                        variogram = out$vgram)

  
  # Semi-isotropic variogram
  # Average the variogram for cells with different directions, but same
  # row/col separation. Since out$ind has only positive column separations,
  # we need to average rows with the same absolute value of rows and the same
  # value of columns
  dat.halfheat <- aggregate(dat.aniso$z,
                            by = abs(dat.aniso[, c('x', 'y')]), mean)
  names(dat.halfheat)[3] <- 'z'

  #   var.half.matrix <- tapply(out$vgram.full, as.data.frame(abs(out$ind)), mean)
  #   data.frame(expand.grid(lapply(dimnames(var.half.matrix),
  #                                 as.numeric)),
  #              z = as.vector(var.half.matrix))
  #   image(var.half.matrix)

  
  # The perspective plot can only be done with base graphics
  # I do the computations here, and pass arguments for the print method
  ind <- lapply(dat.halfheat[, c('x', 'y')], as.factor)
  dat.mhalf <- matrix(NA, nrow = nlevels(ind$x), ncol = nlevels(ind$y))
  dat.mhalf[sapply(ind, as.numeric)] <- dat.halfheat$z
  
  #   jet.colors <- colorRampPalette( c("blue", "green") )
  col.f <- colorRampPalette( c(breedR.getOption('col.seq')[1],
                               breedR.getOption('col.seq')[2]) )
  nbcol = 100
  colours <- col.f(nbcol)
  
  # Computes the facet's means without taking into account NAs
  meanmat.f <- function(x) {
    nx = nrow(x)
    ny = ncol(x)
    tmp <- cbind(as.vector(x[-1, -1]),
                 as.vector(x[-1, -ny]),
                 as.vector(x[-nx, -1]),
                 as.vector(x[-nx, -ny]))
    
    mv <- apply(tmp, 1, mean, na.rm = TRUE)
    return(matrix(mv, nrow = nx - 1))
  }
  
  facetcol <- cut(meanmat.f(dat.mhalf), nbcol)
  
  # I can't store the plot itself in an object
  # but I can stor the arguments for persp() in a list
  dat.halfpersp <- list(x = as.numeric(levels(ind$x)),
                        y = as.numeric(levels(ind$y)),
                        z = dat.mhalf,
                        zlim = c(0, max(dat.mhalf, na.rm = TRUE)),
                        col = colours[facetcol], phi = 30, theta = -35,
                        box = TRUE,
                        axes = TRUE,
                        ticktype = 'detailed',
                        xlab = 'row disp.',
                        ylab = 'col disp.',
                        zlab = '')

  #   # DEBUG
  #   # How colour assignement works
  #   nx <- ny <- 4
  #   tdat <- matrix(NA, nrow = nx, ncol = ny)
  #   tz <- 1:(nx*(ny+1)/2)
  #   tdat[c(1:7, 9:10, 13)] <- tz
  #   tcol <- col.f(100)
  #   # Compute the z-value at the facet centres
  #   # This gives problems with NAs 
  #   tzfacet <- tdat[-1, -1] + tdat[-1, -ny] + tdat[-nx, -1] + tdat[-nx, -ny]
  #   # Better use the values at the vertex closest to the origin
  #   tzfacet <- tdat[-nx, -ny]
  #   # Performing the mean with matrices with NAs
  #   tzfacet <- meanmat.f(tdat)
  #   tfcol <- cut(tzfacet, 50)
  #   persp(x = 1:4, y = 1:4, z = tdat, col = tcol[tfcol], phi = 30, theta = -35)
  
  # Visualization
  
  #   # Residuals
  #   image(dat)
  #   ggplot(data.frame(coord1, z = z), aes(x, y, fill = z)) +
  #     geom_raster() + 
  #     coord_fixed() +
  #     scale_fill_gradient(low='green', high='red')

  
  ans <- list(isotropic = dat.iso,
              perspective = dat.halfpersp,
              heat = dat.halfheat,
              anisotropic = dat.aniso,
              plot = plot)
  class(ans) <- 'breedR.variogram'
  return(ans)
}

#' @method print breedR.variogram
#' @import ggplot2
#' @export
print.breedR.variogram <- function(x, ...) {
  
  # Compute the relevant plots only
  # Except for 'perspective' that connot be precomputed
  if(x$plot == 'isotropic' | x$plot == 'all') {
    p.iso <-   ggplot(x$isotropic,
                      aes(distance, variogram)) +
      geom_point() +
      geom_line() +
      stat_smooth(se = FALSE, method = 'auto')
  }

  if(x$plot == 'heat' | x$plot == 'all') {
    p.heat <-   spatial.plot(x$heat, scale = 'seq')
  }
  
  if(x$plot == 'anisotropic' | x$plot == 'all') {
    p.aniso <-   spatial.plot(x$anisotropic, scale = 'seq')
  }
  
  if(x$plot == 'all') {

    # require(grid)
    if (!requireNamespace("grid", quietly = TRUE)) {
      stop("Package grid needed for plotting all variograms at once. Please install it, or plot them one by one.",
           call. = FALSE)
    }
    
    plot.new()
    # split the graphic window in four and reduce margins
    op <- par(mfrow = c(2, 2), mar = rep(0.5, 4))
    par(mfg = c(2, 1))
    
    # define vewports
    vp.BottomRight <- grid::viewport(height = grid::unit(.5, "npc"),
                                     width = grid::unit(0.5, "npc"), 
                                     just = c("left","top"), 
                                     y = 0.5, x = 0.5)
    vp.TopLeft <- grid::viewport(height = grid::unit(.5, "npc"),
                                 width = grid::unit(0.5, "npc"), 
                                 just = c("right","bottom"), 
                                 y = 0.5, x = 0.5)
    vp.TopRight <- grid::viewport(height = grid::unit(.5, "npc"),
                                  width = grid::unit(0.5, "npc"), 
                                  just=c("left","bottom"), 
                                  y = 0.5, x = 0.5)
    
    # print the four plots
    suppressMessages(print(p.iso, vp=vp.TopLeft))
    print(p.heat, vp=vp.TopRight)
    do.call('persp', args = x$perspective)
    print(p.aniso, vp=vp.BottomRight)

    # Return graphical device to default state
    par(op)
  }
  
  if(x$plot == 'isotropic') suppressMessages(print(p.iso))
  if(x$plot == 'anisotropic') print(p.aniso)
  if(x$plot == 'heat') print(p.heat)
  if(x$plot == 'perspective') do.call('persp', args = x$perspective)
  if(x$plot == 'none') str(x, 1)
}