#' Compare two or more ggplots of the same kind
#' 
#' This function presents several ggplots of the same type side by side
#' under the same scale, while keeping annotations.
#' 
#' The names of the objects in the list will be used for facet labels.
#' 
#' @param plots List of ggplots with meaningful names
#' 
#' @usage compare.plots(plots)
#' @import ggplot2
#' @export compare.plots
compare.plots <- function(plots) {
  # require(plyr)
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop("Package plyr needed for comparing plots. Please install..",
         call. = FALSE)
  }
  
  # Use the same parameters as one of the plots, and add a facet
  # Thus plots need to have "compatible" names
  p <- plots[[1]]
  
  # Aggregate datasets and substitute data
  # http://docs.ggplot2.org/current/gg-add.html
  if( is.null(names(plots)) ) {
    names(plots) <- paste('p', seq_along(plots), sep = '')
  }
  tmpdat <- plyr::ldply(plots, function(x) x$data)
  # Keep the order of the plots
  tmpdat <- transform(tmpdat,
                      .id = factor(.id, levels = unique(tmpdat$.id)))
  p <- p %+% tmpdat
  # Annotations
  extract.text.data <- function(plot) {
    # identify the geom_text layer
    lab.idx <- which(plyr::laply(plot$layers,
                           function(x) x$geom$objname == 'text'))
    if(length(lab.idx) == 0) return(NULL)
    data.frame(layer = lab.idx,
               plot$layers[[lab.idx]]$data,
               parse = plot$layers[[lab.idx]]$geom_params$parse)
  }
  text.data <- plyr::ldply(plots, extract.text.data)

  # If there are annotations ...
  if( nrow(text.data) > 0 ) {
    # Remove the original geom_text layer
    p$layers[[text.data[1, 'layer']]] <- NULL
    p <- p + geom_text(aes_string("x", "y", label = "lab"),
                       data = text.data,
                       parse = any(text.data$parse))    
  }
  
  # Include all annotations and facets
  p <- p + facet_grid("~ .id")
  p
  #   
  #   ggplot(tmpdat, aes(irow, icol)) + 
  #     geom_tile(aes(fill = value)) + 
  #     scale_fill_gradient(low = 'green', high = 'red') +
  #     facet_wrap(~ .id)
  #   qplot(irow, icol, fill = value, data = tmpdat, geom = 'tile') + facet_wrap(~.id)
}

#' Plot an spatially arranged continuous variable
#' 
#' @param dat A 3-column data.frame with names 'x', 'y' and 'z' where the first 
#'   two are the spatial coordinates, and 'z' is the value to be represented
#' @param scale Character. 'divergent' represents positive and negative values 
#'   with different colours. 'sequential' uses a gradient scale of two colours.
#' @import ggplot2
spatial.plot <- function(dat, scale = c('divergent', 'sequential')) {
  
  scale <- match.arg(scale)
  
  dat <- as.data.frame(dat)
  cn <- names(dat)
  
  if( !all(c('x', 'y', 'z') %in% cn) ) {
    ggcl <- paste('ggplot2::ggplot(dat, aes(',cn[1], ',', cn[2], ')) + geom_raster(aes(fill = ', cn[3], '))')
    p <- eval(parse(text = ggcl))
  } else {
    p <- ggplot2::ggplot(dat, aes_string("x" , "y")) +
      geom_raster(aes_string(fill =  "z"))
  }
  
    
  p <- p + coord_fixed()
  
  # Tool to extract hex-codes of colours
  #   scale_colour_brewer(type = 'seq', palette = 'Oranges')$palette(8)
  #   http://blog.ggplot2.org/post/24607351280/choosing-colour-palettes-part-ii-educated-choices
  
  p <- switch(scale,
              divergent = p + scale_fill_gradient2(low  = breedR.getOption('col.div')[1],
                                                   high = breedR.getOption('col.div')[2],
                                                   space = 'Lab'),
              sequential = p + scale_fill_gradient(low  = breedR.getOption('col.seq')[1],
                                                   high = breedR.getOption('col.seq')[2])
  #                 sequential = p + scale_fill_gradient(low = "#FFF7FB", high = "#034E7B")
  #                 sequential = p + scale_fill_gradientn(colours = terrain.colors(10))
  #                 sequential = p
  #                 sequential = p + scale_fill_gradient(low = 'black', high = 'white')
  #                 sequential = p + scale_fill_gradient(low = scales::muted('green'),
  #                                                      high = scales::muted('red'))
              )
  p
}
