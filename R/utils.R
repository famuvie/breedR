## Internal utility functions
## Not exported

#' lme4-style formulas
#' 
#' Transform the separated fixed and random formulas
#' into the single formula with lme4 syntaxis
lme4_fml <- function(fix, rnd, rm_int = TRUE) {
  rnd.terms <- attr(terms(rnd), 'term.labels')
  rnd.terms.lme4 <- paste('(1|', rnd.terms, ')', sep ='')
  int <- ifelse(rm_int, '-1', '')
  rnd.upd <- paste('~ .', int, paste('+', rnd.terms.lme4, collapse = ' '))
  fml.lme4 <- update(fix, rnd.upd)
  return(fml.lme4)
}


breedR.is.element <- function(name, alist) {
  ## return TRUE if element with name NAME is a member of LIST and
  ## the value is non null and not NA.
  if (any(names(alist) == name)) {
    idx = which(names(alist) == name)
    if (!is.null(alist[[idx]]) && !suppressWarnings(is.na(alist[[idx]]))) {
      return (TRUE)
    } else {
      return (FALSE)
    }
  } else {
    return (FALSE)
  }
}

breedR.get.element <-  function(name, alist) {
  if (breedR.is.element(name, alist)) {
    return (alist[[which(names(alist) == name)]])
  } else {
    return (NULL)
  }
}


#' Fit some model
breedR.result <- function(...) {
  res  <- suppressWarnings(remlf90(fixed  = phe_X ~ gg,
                                   genetic = list(model = 'add_animal', 
                                                  pedigree = globulus[,1:3],
                                                  id = 'self'), 
                                   spatial = list(model = 'AR', 
                                                  coord = globulus[, c('x','y')],
                                                  rho = c(.85, .8)), 
                                   data = globulus,
                                   ...))
  return(res)
}

#' Geometric mean
gmean <- function(x) {
  logx <- log(x)
  finite.logx <- is.finite(logx)
  if( !all(finite.logx) ) {
    warning('Removing zeroes for geometric mean')
    logx <- logx[finite.logx]
  }
  
  return(exp(mean(logx)))
}


#' BreedR binaries
#' 
#' Return the path to breedR binaries.
#' Path is different in each platform, but not architecture.
`breedR.bin.builtin` = function()
{
  #   if (breedR.os("mac")) {
  #     fnm = system.file(paste("bin/mac/", breedR.os.32or64bit(), "bit/breedR", sep=""), package="breedR")
  #   } else if (breedR.os("linux")) {
  #     fnm = system.file(paste("bin/linux/breedR", breedR.os.32or64bit(), sep=""), package="breedR")
  #   } else if (breedR.os("windows")) {
  #     fnm = system.file(paste("bin/windows/breedR", breedR.os.32or64bit(), ".exe", sep=""), package="breedR")
  #   } else {
  #     stop("Unknown OS")
  #   }
  
  if( breedR.os.type() != 'else') {
    fnm <- system.file('bin', breedR.os.type(), package='breedR')
  } else {
    stop("Unknown platform")
  }
  
  if (file.exists(fnm)) {
    return (fnm)
  } else {
    stop(paste("breedR installation error; no such file", fnm))
  }
}


#' Determine the user's home directory
#' 
#' Relies on \code{Sys.getenv('HOME')}, or under windows, on
#' \code{Sys.getenv("USERPROFILE"))} changing backslashes to slashes.
`breedR.get.HOME` = function()
{
  return (as.character(ifelse(breedR.os("windows"),
                              gsub("\\\\", "/", Sys.getenv("USERPROFILE")),
                              Sys.getenv("HOME"))))
}

#' Determine the user name
`breedR.get.USER` = function()
{
  u = ""
  for (U in c("USER", "USERNAME", "LOGNAME")) {
    u = Sys.getenv(U)
    if (u != "")
      break;
  }
  if (u == "")
    u = "UnknownUserName"
  
  return (as.character(u))
}
