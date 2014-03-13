#' Compare two or more ggplots of the same kind
#' 
#' This function presents several ggplots of the same type side to side
#' under the same scale, while keeping annotations.
#' 
#' @param plots List of ggplots with meaningful names
#' 
#' The names of the objects in the list will be used for facet labels.
#' @export
compare.plots <- function(plots) {
  require(plyr)
  # Use the same parameters as one of the plots, and add a facet
  # Thus plots need to have "compatible" names
  p <- plots[[1]]
  
  # Aggregate datasets and substitute data
  # http://docs.ggplot2.org/current/gg-add.html
  tmpdat <- ldply(plots, function(x) x$data)
  # Keep the order of the plots
  tmpdat <- transform(tmpdat,
                      .id = factor(.id, levels = unique(tmpdat$.id)))
  p <- p %+% tmpdat
  # Annotations
  extract.text.data <- function(plot) {
    # identify the geom_text layer
    lab.idx <- which(laply(plot$layers,
                           function(x) x$geom$objname == 'text'))
    if(length(lab.idx) == 0) return(NULL)
    data.frame(layer = lab.idx,
               plot$layers[[lab.idx]]$data,
               parse = plot$layers[[lab.idx]]$geom_params$parse)
  }
  text.data <- ldply(plots, extract.text.data)
  
  # If there are annotations ...
  if(nrow(text.data) > 0) {
    # Remove the original geom_text layer
    p$layers[[text.data[1, 'layer']]] <- NULL
    p <- p + geom_text(aes(x, y, label = lab),
                       data = text.data,
                       parse = any(text.data$parse))    
  }
  
  # Include all annotations and facets
  p <- p + facet_grid(~ .id)
  p
  #   
  #   ggplot(tmpdat, aes(irow, icol)) + 
  #     geom_tile(aes(fill = value)) + 
  #     scale_fill_gradient(low = 'green', high = 'red') +
  #     facet_wrap(~ .id)
  #   qplot(irow, icol, fill = value, data = tmpdat, geom = 'tile') + facet_wrap(~.id)
}

