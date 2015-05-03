#' @export
all.equal.remlf90 <- function(target, current, check.attributes = TRUE, use.names= TRUE, ...) {
  
  ## Exclude some elements of the remlf90 objects
  ## that do not include numerical results
  ## and are not usually interesting for comparison
  exclusions <- c('call', 'method', 'components', 'effects', 'mf', 'mm', 'reml')
  thin_tg <- target[-match(exclusions, names(target))]
  thin_cr <- current[-match(exclusions, names(target))]
  
  all.equal.list(thin_tg, thin_cr, check.attributes = check.attributes, use.names = use.names, ...)
}