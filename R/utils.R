## Internal utility functions
## Not exported

# lme4-style formulas
# 
# Transform the separated fixed and random formulas
# into the single formula with lme4 syntaxis
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


# Fit some model
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

# Geometric mean
gmean <- function(x) {
  logx <- log(x)
  finite.logx <- is.finite(logx)
  if( !all(finite.logx) ) {
    warning('Removing zeroes for geometric mean')
    logx <- logx[finite.logx]
  }
  
  return(exp(mean(logx)))
}


# BreedR binaries
# 
# Return the default path to breedR binaries checking its existence.
breedR.bin.builtin  <- function()
{

  if( breedR.os.type() != 'else') {
    bindir <- 'bin'
    if (breedR.os('windows'))
      bindir <- file.path(bindir, paste0(breedR.os.32or64bit(), 'bit'))
    fnm <- system.file(bindir, package='breedR')
  } else {
    stop("Unknown platform")
  }
  
  if (file.exists(fnm)) {
    return (fnm)
  } else {
    stop(paste("breedR installation error; no such directory", fnm))
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

# Convert a incidence matrix specified in 8+8 columns format
# the first 8 are coefficients, and the last 8 are columns
# into a sparse matrix format
matrix.short16 <- function(M) {
  coef = M[, 1:8]
  neig = M[, 8+1:8]
  
  n <- nrow(coef)
  p <- max(neig, na.rm = TRUE)
  
  i <- rep(1:n, 8)
  j <- as.vector(neig)
  x <- as.vector(coef)
  
  rm.idx <- which(x==0)
  stopifnot( all(x[rm.idx] == 0) )
  
  Z <- Matrix::spMatrix(nrow = n, ncol = p,
                        i = i[-rm.idx],
                        j = j[-rm.idx],
                        x = x[-rm.idx])
  return(Z)
}

#' 'Splat' arguments to a function
#' 
#' Borrowed from \code{\link[plyr]{splat}}
splat <- function (flat) {
  function(args, ...) {
    do.call(flat, c(args, list(...)))
  }
}

# given a list of data.frames, extract a given column
# from each of the data.frames into a matrix.
# Optionally drop into a vector if dimension = 1.
# Used in ranef.remlf90 and fixef.remlf90 to extract
# trait-wise predictions of effects
ldf2matrix <- function(x, vname, drop = TRUE) {
  ## All dataframes (ntraits) are of the same size (nlevels x (value, s.e.))
  ## ensure a matrix, even if nlevels = 1
  ans <- do.call(cbind, lapply(x, `[[`, vname))
  rownames(ans) <- rownames(x[[1]])
  if (drop) ans <- drop(ans)
  return(ans)
}


# Extract values and standard errors from lists
# of effects estimates (or predictions)
# x is a list of effects, where each element is a trait-wise
# list of data.frames with columns 'value' and 's.e.'
get_estimates <- function(x) {
  values <- lapply(x, ldf2matrix, 'value')
  se <- lapply(x, ldf2matrix, 's.e.')
  ans <- mapply(function(gvl, gse) structure(gvl, se = gse),
                values, se, SIMPLIFY = FALSE)
  return(ans)
}

# combine sub-effect names and trait names
# trait names within sub-effect names
# bl1 bl2 bl3 + y1 y2 = bl1.y1 bl1.y2 bl2.y1 bl2.y2 bl3.y1 bl3.y2
# names_effect(paste0("bl", 1:3), paste0("y", 1:2))
# "bl1.y1" "bl1.y2" "bl2.y1" "bl2.y2" "bl3.y1" "bl3.y2"
# names_effect(paste0("bl", 1:3), NULL)
# "bl1" "bl2" "bl3"
# names_effect(NULL, paste0("y", 1:2))
# "y1" "y2"
# names_effect(NULL, NULL)
# NULL
names_effect <- function(inner = NULL, outer = NULL) {

  ans <- inner
  if (length(outer) > 1) {
    if (!is.null(inner)) {
      ans <- apply(
        expand.grid(
          outer, inner,
          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
        )[2:1],
        1, paste, collapse = "."
      )
    } else {
      ans <- outer
    }
  }
  return(ans)
}


## component names
vcnames <- function(efname, efdim, trnames = trait_names) {
  
  dim_subtrait <- efdim/ifelse(is.null(trnames), 1, length(trnames))
  if (dim_subtrait > 1) 
    efname <- paste(efname, seq_len(dim_subtrait), sep = "_")
  diag_names <- names_effect(efname, trnames)
  
  ## matrix components include variances and covariances
  ans <- outer(diag_names, diag_names, paste, sep = "_")
  diag(ans) <- diag_names
  return(ans[upper.tri(ans, diag = TRUE)])
}
