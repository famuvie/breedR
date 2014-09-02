## Export: breedR.setOption breedR.getOption

#' Set and get global options for breedR
#'
#' Set and get global options for breedR. The options are stored in the variable 
#'\code{breedR.options} in the \code{.GlobalEnv}-environment, and will therefore persist during the session. If you want to set some options permanently do it in a file names \code{.breedRrc} in your home directory. See Examples.
#'
#' @param ... Option and value,  like \code{option=value} or \code{option, 
#'  value}; see the Examples.
#' @param option The option to get. If missing or \code{NULL}, then 
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
#'  \code{cygwin.home}: the user's home in the Cygwin installation [Remote 
#'  computing for Windows only]
#'  
#'  \code{ssh.auth.sock}: the ssh bind-adress (value of $SSH_AUTH_SOCK in the 
#'  Cygwin-shell). [Remote computing for Windows only]
#'  
#'  \code{remote.host}: (IP or DNS) address of a Linux server for remote 
#'  computing
#'  
#'  \code{remote.user}: user name with ssh-keys access to the server (see 
#'  details)
#'  
#'  \code{remote.port}: port for ssh connection
#'  
#'  \code{remote.bin}: path to the binaries directory in the remote installation
#'  of breedR. Usually the output of \code{system.file('bin/linux',
#'  package='breedR')} in a remote server's R-session.
#'  
#'  \code{ssh.options}: ssh options. You shouldn't need to change this.
#'  
#' @details
#'
#' Sequential scales are used for variables not necessarily centered such as a 
#'response variable, or the fitted values of a model. The colour scale is built 
#'as a gradient between two extreme colours which are specified as hex codes or 
#'colour names in the option \code{col.seq}.
#'
#' Diverging scales are used for plots such as residuals, centered (hopefully) 
#'around zero, with positive and negative values represented with different 
#'colours whose intensity is linked to the magnitude. The option \code{col.div} 
#'is a vector of two hex codes or colour names of the most intense colours.
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
#' 
#' \dontrun{
#' # Set up some options permanently by writing the following in $HOME/.breedRc
#' remote.host = '123.45.678.999'
#' remote.user = 'uname'
#' remote.bin  = 'remote/path/to/breedR/bin/linux'
#' }
#' @export breedR.setOption breedR.getOption



breedR.getOption <- function(option = c("ar.eval",
                                        "breedR.bin",
                                        "splines.nok",
                                        "default.initial.variance",
                                        "col.seq",
                                        "col.div",
                                        "cygwin",
                                        "cygwin.home",
                                        "ssh.auth.sock",
                                        "remote.host",
                                        "remote.user",
                                        "remote.port",
                                        "remote.bin",
                                        "ssh.options")) {
  envir <- breedR.get.breedREnv()

  # opt: options set by the user
  option <- match.arg(option, several.ok = TRUE)
  if( exists("breedR.options", envir = envir) ) {
    opt = get("breedR.options", envir = envir)
  } else {
    opt = list()
  }

  # default options in the package
  default.opt = list(
    breedR.bin  = breedR.bin.builtin(),
    ar.eval     = c(-8, -2, 2, 8)/10,
    splines.nok = quote(determine.n.knots),
    default.initial.variance = 1,
    col.seq = c('#034E7B', '#FDAE6B'),
    col.div = c('#3A3A98FF', '#832424FF'),
    cygwin = 'C:/cygwin',
    cygwin.home = paste("/home/", breedR.get.USER(), sep=""),
    ssh.auth.sock = paste("/tmp/ssh-auth-sock-", breedR.get.USER(), sep=""),
    remote.host = "your.computing.server",
    remote.user = "yourusername",
    remote.port = 22,
    remote.bin  = "path_to/breedR/bin/linux",
    ssh.options = "-x -o BatchMode=yes -o TCPKeepAlive=yes -e none"
  )
  
  # overriden default through a configuration file
  if( file.exists(file.path(breedR.get.HOME(), '.breedRrc')) ) {
    rc <- new.env()
    sys.source(file.path(breedR.get.HOME(), '.breedRrc'), rc)
    rc <- as.list(rc)
    
    # Check that all options in .breedRrc are legal options
    if( !all(names(rc) %in% names(default.opt)) ) {
      unknown_opt <- ls(rc)[which(!ls(rc) %in% names(default.opt))]
      stop(paste('Option(s) not recognized in ~/.breedRrc:',
                 paste(unknown_opt, collapse = ', ')))
    }
    
    # Override defaults from .breedRrc
    default.opt[names(rc)] <- rc
  }
  
  # if asked for no particular option, then show all
  if (missing(option) | is.null(option))
    option <- names(default.opt)
  #     return(default.opt)   
  #     stop("argument is required.")
  
  # Return current values for the requested options if set,
  # or default values if not set.
  res = vector('list', length(option))
  names(res) <- option
  for (i in seq_along(option)) {
    if (breedR.is.element(option[i], opt)) {
      res[[option[i]]] <- breedR.get.element(option[i], opt)
    } else {
      res[[option[[i]]]] <- breedR.get.element(option[i], default.opt)
    }
  }
  
  # Return only the value if only one option was requested
  # (instead of a list of length one)
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
                                                "ssh.auth.sock",
                                                "remote.host",
                                                "remote.user",
                                                "remote.port",
                                                "remote.bin",
                                                "ssh.options"),
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
