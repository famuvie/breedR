
#' Empirical variograms of residuals
#' 
#' Computes isotropic and anisotropic empirical variograms from the residuals of
#' a breedR model.
#' 
#' This assumes that there is at most one observation per spatial location. 
#' Otherwise, are observations measured at different times? should a
#' spatial-temporal variogram be fitted?
#' 
#' @importFrom fields vgram.matrix
#' @export
variogram <- function(x, R, plot = c('all', 'isotropic', 'anisotropic', 'perspective', 'heat', 'none')) {
  
  if( !inherits(x, 'breedR') ) 
    stop('This function works on models fitted with breedR')
  
  plot <- match.arg(plot)
  # We need to place the residuals in matrix form
  # where rows and columns have the same spacing
  coord <- coordinates(x)   # Coordinates in distance units (e.g. m)
  resid <- residuals(x)
  
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
  dat <- sparseMatrix(i = as.integer(coord1[[1]]),
                      j = as.integer(coord1[[2]]),
                      x = resid
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
  col.f <- colorRampPalette( c(breedR.getOption('colours.seq')[1],
                               breedR.getOption('colours.seq')[2]) )
  nbcol = 100
  colours <- col.f(nbcol)
  facetcol <- cut(dat.halfheat$z, nbcol)
  
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

  # Visualization
  
  #   # Residuals
  #   image(dat)
  #   ggplot(data.frame(coord1, z = resid), aes(x, y, fill = z)) +
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

#' @S3method print breedR.variogram
#' @export
print.breedR.variogram <- function(x) {
  
  # Compute the relevant plots only
  # Except for 'perspective' that connot be precomputed
  if(x$plot == 'isotropic' | x$plot == 'all') {
    p.iso <-   ggplot(x$isotropic,
                      aes(distance, variogram)) +
      geom_point() +
      stat_smooth(se = FALSE, method = 'auto')
  }

  if(x$plot == 'heat' | x$plot == 'all') {
    p.heat <-   spatial.plot(x$heat, scale = 'seq')
  }
  
  if(x$plot == 'anisotropic' | x$plot == 'all') {
    p.aniso <-   spatial.plot(x$anisotropic, scale = 'seq')
  }
  
  if(x$plot == 'all') {
    require(grid)
    
    plot.new()
    # split the graphic window in four
    op <- par(mfrow = c(2, 2), mar = rep(0.5, 4))
    par(mfg = c(2, 1))
    
    # define vewports
    vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                               just=c("left","top"), 
                               y=0.5, x=0.5)
    vp.TopLeft <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("right","bottom"), 
                           y=0.5, x=0.5)
    vp.TopRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                            just=c("left","bottom"), 
                            y=0.5, x=0.5)
    
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