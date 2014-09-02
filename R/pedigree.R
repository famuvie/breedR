#' Build pedigree
#' 
#' Builds a full pedigree out of observed data with sorting and recoding.
#' 
#' A full pedigree requires that all the individual codes for sires or dams are 
#' present as individuals themselves, possibly with unknown parents. Therefore, 
#' for using it in a statistical model, it is necessary to \emph{complete} the 
#' pedigree by introducing new individuals with unknown parents.
#' 
#' Furthermore, the codes must be sorted in ascending and consecutive order 
#' begining from 1, and the offspring must follow parents. All this is checked, 
#' and the pedigree is reordered and recoded if needed.
#' 
#' If recoding is needed, the function issues a warning and an attribute 'map' 
#' is attached to the pedigree, such that \code{map[i] = j} means that code 
#' \code{i} was renumbered as \code{j}. Therefore, if \code{x} is a vector with 
#' original codes, \code{map[x]} gives the new codes. Conversely, \code{match(y,
#' map)} back-transforms to original codes.
#' 
#' @param x if given, a vector of length 3 with indices or names of columns in 
#'   \code{data} corresponding to \code{self}, \code{sire} and \code{dam} codes
#' @param self index or column name in \code{data} with codes of sires
#' @param sire index or column name in \code{data} with codes of individuals
#' @param dam index or column name in \code{data} with codes of dams
#' @param data a dataframe or a list to take the individual codes from
#' @return A well-formed 'pedigree'-class object. Possibly sorted and recoded.
#' @seealso \code{link{check_pedigree}}
#' @examples
#' # Founders are missing in the globulus dataset
#' data(globulus)
#' check_pedigree(globulus[,c('self', 'dad', 'mum')])
#' # build_pedigree completes the missing information
#' ped <- build_pedigree(c('self', 'dad', 'mum'), data = globulus)
#' check_pedigree(ped)
#' @export
build_pedigree <- function(x, self = x[[1]], sire = x[[2]], dam = x[[3]], data) {
  ### Checks
  missing_terms <- c(self = missing(self), 
                     sire = missing(sire), 
                     dam = missing(dam))
  if( missing(x) & any(missing_terms) )
    stop("Provide either 'x' or the triplet 'self', 'sire' and 'dam'")
  if( missing(data) ) stop("missing 'data'" )
  if( !missing(x) ) {
    if( length(x) != 3 ) stop("'x' must be of length 3" )
    if( !any(missing_terms ) )
      warning(paste(paste(names(which(missing_terms)), collapse = ' and '), 
                    "will override 'x'"))
  }
  if( !all(sapply(data[c(self, sire, dam)], is.numeric)) )
    stop("The pedigree codes must be numeric")
  
  # Extract the relevant columns of the dataframe
  ped <- as.data.frame(data)[c(self, sire, dam)]
  names(ped) <- c('self', 'sire', 'dam')

  # 0 is to be interpreted as unknown parent, not as an individual code
  # recode it as NA
  ped$sire[which(ped$sire == 0)] <- NA
  ped$dam[which(ped$dam == 0)] <- NA
  
  # Codes of sires or dams not present as self
  missing_codes <- unique(c(ped$sire[! ped$sire %in% ped$self],
                            ped$dam[! ped$dam %in% ped$self]))

  # Except NA which is not an identified individual
  missing_codes <- missing_codes[missing_codes>0 & !is.na(missing_codes)]
  
  n.add <- length(missing_codes)
  
  # Complete the pedigree with missing codes
  pedx <- rbind(data.frame(self = sort(missing_codes),
                           sire = rep(NA, n.add),
                           dam  = rep(NA, n.add)),
                ped)

  # Order by identity
  pedx <- pedx[order(pedx$self), ]
  ### Check pedigree rules
  checks <- check_pedigree(pedx)
  
  # The pedigree should now be full and sorted
  stopifnot(all(checks[c('full_ped', 'codes_sorted')]))
  
  # If any other rule is not met, reorder and recode
  # Build a map from position to code
  if( !all(checks) ) {
    ord <- pedigree::orderPed(pedx)
    map <- rep(NA, max(pedx[, 1]))
    map[pedx[order(ord), 1]] <- 1:nrow(pedx)
    pedx <- as.data.frame(sapply(pedx, function(x) map[x]))
    pedx <- pedx[order(pedx$self), ]
    warning("The pedigree has been recoded. Check attr(ped, 'map').")
  }
  
  # pedx should pass all checks
  stopifnot(all(check_pedigree(pedx)))

  out <- pedigreemm::pedigree(sire = pedx$sire, dam = pedx$dam, label = pedx$self)
  if( !all(checks)) attr(out, 'map') <- map
  return(out)
}

#' Check the rules for a well-formed pedigree
#' 
#' This function performs sanity checks on a pedigree.
#' 
#' Rules are: the pedigree must be full (i.e. all the codes in the pedigree must
#' appear once in the first column); the codes of the offspring must follow 
#' their parent's codes; the identity codes must be sorted in ascending order
#' and this order must be consecutive begining from 1.
#' 
#' If any check fails, the pedigree must be recoded/reordered to be used in 
#' analysis. \code{build_pedigree(1:3, data = ped)} should fix it.
#' 
#' @param ped a 'pedigree'-class object, a dataframe, a matrix, or anything that
#'   can be coerced into a 3-column, numeric dataframe.
#'   
#' @return a named logical vector with the results of the checks
#' @seealso \code{\link{build_pedigree}}
#' @examples
#' # A well-formed pedigree
#' ped_ok <- data.frame(id   = 1:6,
#'                      dam  = c(NA,NA,1,1,4,4),
#'                      sire = c(NA,NA,2,2,3,5))
#' check_pedigree(ped_ok)  # passes all checks
#' 
#' # Sometimes founders are missing
#' (ped_notfull <- ped_ok[-c(1:2),])
#' check_pedigree(ped_notfull)  # fails full_ped
#' 
#' # Sometimes codes of parents are greater than their offspring
#' sw_23 <- c(1, 3, 2, 4:6)
#' (ped_messed <- as.data.frame(sapply(ped_ok, function(x) sw_23[x]))[sw_23,])
#' check_pedigree(ped_messed)        # fails offsp_follows
#' 
#' # Sometimes the pedigree is unordered
#' (ped_unordered <- ped_ok[sw_23, ])
#' check_pedigree(ped_unordered)   # fails codes_sorted
#' 
#' # Finally, sometimes codes are just not consecutive
#' (ped_notconsec <- transform(ped_ok, id = c(1:5, 10)))
#' check_pedigree(ped_notconsec)  # fails codes_consecutive
#' 
#' # Everything is fixed with build_pedigree()
#' check_pedigree(build_pedigree(1:3, data = ped_notfull))
#' check_pedigree(build_pedigree(1:3, data = ped_messed))
#' check_pedigree(build_pedigree(1:3, data = ped_unordered))
#' check_pedigree(build_pedigree(1:3, data = ped_notconsec))
#' @export
check_pedigree <- function(ped) {
  
  # ped might be a matrix, a pedigreemm::pedigree
  # or anything that can be coerced to a dataframe
  ped <- as.data.frame(ped)
  N   <- nrow(ped)

  # but at the end of the day, it has to have 3 columns with numbers
  stopifnot(ncol(ped) == 3L)
  stopifnot(all(sapply(ped, is.numeric)))
  
  # The first column is assumed to be the individuals' codes
  # and the other two, the codes of both progenitors
  
  # Parent codes should appear in the first column
  # This removes 0 and NA's as codes
  parent_codes <- sort(unique(stack(ped[,2:3])$values[which(ped[,2:3] > 0)]))
  full_ped <- all(parent_codes %in% ped[, 1])
  
  # Offspring follow parents (in codes)
  offsp_follows <- all(apply(ped, 1, function(x) all(x[1] > x[2:3], na.rm = TRUE)))
  
  # Codes are sorted
  ord_codes <- order(ped[, 1])
  codes_sorted  <- identical(ord_codes, 1:N)
  
  # Codes are consecutive
  rg_codes <- range(ped[, 1])
  codes_consec <- identical(as.integer(ped[ord_codes, 1]), 
                            seq(rg_codes[1], rg_codes[2]))
  
  return(c(full_ped     = full_ped,
           offsp_follows = offsp_follows,
           codes_sorted = codes_sorted,
           codes_consecutive = codes_consec))
}

# # Tests for sorting as graphs
# # A pedigree is a Directed Acyclic Graph (DAG)
# # and as such, it has a topological order, 
# # where ancestors preceedes descendants
# # There are linear algorithms for finding a topological order.
# # Implementations in R:
# # (CRAN) gRbase::topoSort (calls a C function)
# # (BioC) RBGL::sort interfaces the C++ BGL library of graph algorithms 
# # (CRAN) pedigree::orderPed (calls a C function, doesn't need a graph)
# 
# nV <- 8
# testV <- letters[1:nV]
# edL <- list(a = list(edges = 'd'),
#             b = list(edges = c('d', 'e', 'f')),
#             c = list(edges = 'e'),
#             d = list(edges = 'g'),
#             e = list(edges = 'g'),
#             f = list(edges = NULL),
#             g = list(edges = 'h'),
#             h = list(edges = NULL))
# testg <- graphNEL(nodes = testV, edgeL = edL, edgemode = 'directed')
# plot(testg)
# gRbase::topoSort(testg)
# RBGL::tsort(testg)
# 
# testAM <- as(as(testg, 'matrix'), 'dgTMatrix')
# perm <- sample(nV, nV)
# twistAM <- testAM[perm, perm]
# plot(as(as.matrix(twistAM), 'graphAM'))
# testg_df <- as.data.frame(cbind(1:nV, t(sapply(1:nV, function(x) c((testAM@i + 1)[testAM@j + 1 == x], NA, NA)[1:2]))))
# twist_df <- as.data.frame(cbind(1:nV, t(sapply(1:nV, function(x) c((twistAM@i + 1)[twistAM@j + 1 == x], NA, NA)[1:2]))))
# perm_rev <- pedigree::orderPed(twist_df)
# perm[order(perm_rev)] # Not the identity, but a valid renumbering:
# plot(as(as.matrix(twistAM[order(perm_rev), order(perm_rev)]), 'graphAM'))
# 
# # Better way of twisting and untwisting without using matrices
# # Furthermore, assigning consecutive numbers from 1 to nV
# twist2 <- sapply(testg_df, function(x) order(perm)[x])
# twist2[order(twist2[,1]), ]
