## Export: breedR.setOption breedR.getOption

#' Set and get global options for breedR
#' 
#' Set and get global options for breedR. The options are stored in the variable
#' \code{breedR.options} in the \code{.GlobalEnv}-environment.
#' 
#' @param ... Option and value,  like \code{option=value} or \code{option, 
#'   value}; see the Examples.
#' @param option The option to get. If missing or \code{NULL}, then 
#'   \code{breedR.getOption} will display the current defaults, otherwise, 
#'   \code{option} must be one of
#'   
#'   \code{ar.eval}: numeric vector of values in (-1, 1) where the 
#'   autoregressive parameters should be evaluated if not otherwise specified
#'   
#'   \code{splines.nok}: a function of the number of individuals in a row which 
#'   gives the number of knots (nok) to be used for a splines model, if not 
#'   otherwise specified
#'   
#'   \code{default.initial.variance}: a default value for all variance
#'   components
#'   
#' @name breedR.option
#' @aliases breedR.options breedR.setOption breedR.getOption
#' @examples
#' ## Set default values for the autoregressive parameters
#' breedR.setOption("ar.eval", 3*(-3:3)/10)
#' ## alternative format
#' breedR.setOption(ar.eval = 3*(-3:3)/10) 
#' ## check it 
#' breedR.getOption("ar.eval")
#' @export breedR.setOption breedR.getOption


breedR.getOption <- function(option = c("ar.eval",
                                        "splines.nok",
                                        "default.initial.variance")) {
  default.opt = list(
    ar.eval     = c(-8, -2, 2, 8)/10,
    splines.nok = quote(breedR:::determine.n.knots),
    default.initial.variance = 1
  )
  
  if (missing(option) | is.null(option))
    return(default.opt)
#     stop("argument is required.")
  
  envir = breedR.get.breedREnv()
  
  option = match.arg(option, several.ok = TRUE)
  if (exists("breedR.options", envir = envir))
    opt = get("breedR.options", envir = envir)
  else
    opt = list()
  
  res = c()
  for (i in 1:length(option)) {
    if (breedR.is.element(option[i], opt)) {
      res = c(res, breedR.get.element(option[i], opt))
    } else {
      res = c(res, breedR.get.element(option[i], default.opt))
    }
  }
  
  return (res)
}



## supports the following formats:
##     breedR.setOption("ar.eval", .5)
##     breedR.setOption(ar.eval = .5)
##     breedR.setOption(ar.eval = .5, default.initial.variance = 10)
#' @rdname breedR.option
breedR.setOption <- function(...) {
  
  breedR.setOption.core <-  function(option = c("ar.eval",
                                                "splines.nok",
                                                "default.initial.variance"),
                                     value) {
    envir = breedR.get.breedREnv()
    
    option = match.arg(option, several.ok = FALSE)
    if (!exists("breedR.options", envir = envir))
      assign("breedR.options", list(), envir = envir)
    if (is.character(value)) {
      eval(parse(text = paste("breedR.options$", option, "=", shQuote(value), sep="")),
           envir = envir)
    } else {
      eval(parse(text = paste("breedR.options$", option, "=", ifelse(is.null(value), "NULL", value), sep="")),
           envir = envir)
    }
    return (invisible())
  }
  
  called = list(...)
  len = length(names(called))
  if (len > 0L) {
    for(i in 1L:len) {
      do.call(breedR.setOption.core, args = list(names(called)[i], called[[i]]))
    }
  } else {
    breedR.setOption.core(...)
  }
  return (invisible())
}
