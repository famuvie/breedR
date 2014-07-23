## Export: breedR.setOption breedR.getOption

#'Set and get global options for breedR
#'
#'Set and get global options for breedR. The options are stored in the variable 
#'\code{breedR.options} in the \code{.GlobalEnv}-environment.
#'
#'@param ... Option and value,  like \code{option=value} or \code{option, 
#'  value}; see the Examples.
#'@param option The option to get. If missing or \code{NULL}, then 
#'  \code{breedR.getOption} will display the current defaults, otherwise, 
#'  \code{option} must be one of
#'  
#'  \code{ar.eval}: numeric vector of values in (-1, 1) where the autoregressive
#'  parameters should be evaluated if not otherwise specified
#'  
#'  \code{breedR.bin} path to the breedR binaries
#'  
#'  \code{splines.nok}: a function of the number of individuals in a row which 
#'  gives the number of knots (nok) to be used for a splines model, if not 
#'  otherwise specified
#'  
#'  \code{default.initial.variance}: a default value for all variance components
#'  
#'  \code{col.seq}: a vector with the specification of default extreme breedR
#'  col for sequential scales in spatial quantitative plots. See Details.
#'  
#'  \code{col.div}: a vector with the specification of default extreme breedR
#'  col for diverging scales in spatial quantitative plots. See Details. 
#'  
#'  \code{cygwin}: the home of the Cygwin installation (default "C:/cygwin")
#'  [Remote computing for Windows only]
#'  
#'  \code{cygwin}: the user's home in the Cygwin installation
#'  [Remote computing for Windows only]
#'  
#'  \code{ssh.auth.sock}: the ssh bind-adress (value of $SSH_AUTH_SOCK in the
#'  Cygwin-shell). [Remote computing for Windows only]
#'  
#'@details
#'
#'Sequential scales are used for variables not necessarily centered such as a 
#'response variable, or the fitted values of a model. The colour scale is built 
#'as a gradient between two extreme colours which are specified as hex codes or 
#'colour names in the option \code{col.seq}.
#'
#'Diverging scales are used for plots such as residuals, centered (hopefully) 
#'around zero, with positive and negative values represented with different 
#'colours whose intensity is linked to the magnitude. The option \code{col.div}
#'is a vector of two hex codes or colour names of the most intense colours.
#'
#'@name breedR.option
#'@aliases breedR.options breedR.setOption breedR.getOption
#' @examples
#' ## Set default values for the autoregressive parameters
#' breedR.setOption("ar.eval", 3*(-3:3)/10)
#' ## alternative format
#' breedR.setOption(ar.eval = 3*(-3:3)/10) 
#' ## check it 
#' breedR.getOption("ar.eval")
#'@export breedR.setOption breedR.getOption


breedR.getOption <- function(option = c("ar.eval",
                                        "breedR.bin",
                                        "splines.nok",
                                        "default.initial.variance",
                                        "col.seq",
                                        "col.div",
                                        "cygwin",
                                        "cygwin.home",
                                        "ssh.auth.sock")) {
  envir <- breedR.get.breedREnv()

  option <- match.arg(option, several.ok = TRUE)
  if( exists("breedR.options", envir = envir) ) {
    opt = get("breedR.options", envir = envir)
  } else {
    opt = list()
  }
  
  if ( is.null(opt$breedR.bin) ) {
    breedR.bin = breedR.bin.builtin()
  } else if ( tolower(opt$breedR.bin) == "remote" || tolower(opt$breedR.bin) == "breedR.remote" ) {
    breedR.bin = gsub("\\\\", "/", system.file("bin/remote/breedR.remote", package="breedR"))
  } else
    breedR.bin = opt$breedR.bin
  
  default.opt = list(
    breedR.bin  = breedR.bin,
    ar.eval     = c(-8, -2, 2, 8)/10,
    splines.nok = quote(breedR:::determine.n.knots),
    default.initial.variance = 1,
    col.seq = c('#034E7B', '#FDAE6B'),
    col.div = c('#3A3A98FF', '#832424FF'),
    cygwin = 'C:/cygwin',
    cygwin.home = paste("/home/", breedR.get.USER(), sep=""),
    ssh.auth.sock = paste("/tmp/ssh-auth-sock-", breedR.get.USER(), sep="")
  )
  
  if (missing(option) | is.null(option))
    option <- names(default.opt)
  #     return(default.opt)   
  #     stop("argument is required.")
  
  
  res = list()
  for (i in 1:length(option)) {
    if (breedR.is.element(option[i], opt)) {
      res[[option[i]]] <- breedR.get.element(option[i], opt)
    } else {
      res[[option[[i]]]] <- breedR.get.element(option[i], default.opt)
    }
  }
  
    if(length(res) == 1L) res <- res[[1L]]
  
  return (res)
}



## supports the following formats:
##     breedR.setOption("ar.eval", .5)
##     breedR.setOption(ar.eval = .5)
##     breedR.setOption(ar.eval = .5, default.initial.variance = 10)
#' @rdname breedR.option
breedR.setOption <- function(...) {
  
  breedR.setOption.core <-  function(option = c("ar.eval",
                                                "breedR.bin",
                                                "splines.nok",
                                                "default.initial.variance",
                                                "col.seq",
                                                "col.div",
                                                "cygwin",
                                                "cygwin.home",
                                                "ssh.auth.sock"),
                                     value) {
    
    if(is.list(option)) return(do.call('breedR.setOption', args = option))
    
    envir = breedR.get.breedREnv()
    
    option = match.arg(option, several.ok = FALSE)
    if (!exists("breedR.options", envir = envir))
      assign("breedR.options", list(), envir = envir)
    if (is.character(value)) {
      #       eval(parse(text = paste("breedR.options$", option, "=", shQuote(value), sep="")),
      #            envir = envir)
      envir$breedR.options[[option]] = value
    } else {
      eval(parse(text = paste("breedR.options$", option, "=", ifelse(is.null(value), "NULL", value), sep="")),
           envir = envir)
    }
    return (invisible())
  }
  called = list(...)
  
  # Current options
  op <- lapply(names(called), function(x) do.call('breedR.getOption', args = list(x)))
  names(op) <- names(called)

  len = length(names(called))
  if (len > 0L) {
    for(i in 1L:len) {
      do.call(breedR.setOption.core, args = list(names(called)[i], called[[i]]))
    }
  } else {
    breedR.setOption.core(...)
  }
  return (invisible(op))
}
