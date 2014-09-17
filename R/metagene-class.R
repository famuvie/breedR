#%%%%%%%%%%%%%%%%%%%%#
#%% metagene class %%#
#%% Facundo Muñoz  %%#
#%%%%%%%%%%%%%%%%%%%%#


#' Metagene Data Input
#' 
#' Read the output file (the one with the pedigree. Usually: insim.002) of the 
#' Metagene program, and build an object of class metagene.
#' @param fname file name of the second metagene output file (usually: 
#'   insim.002)
#' @return read.metagene returns the data in the file as an object of class 
#'   \code{metagene}
#' @references \url{http://www.igv.fi.cnr.it/noveltree/}
#' @family metagene
#' @export
read.metagene <- function(fname) {
  file.path = fname
  file      = readLines(con = file.path)
  data_start.idx <- grep('pedigree', file) + 1
  Data <- read.table(file      = file.path,
                     quote     = "",
                     skip      = data_start.idx,
                     col.names = scan(file=file.path,
                                    what=character(),
                                    skip=data_start.idx - 1,
                                    nlines=1,
                                    sep=','))
  t.max = as.integer(strsplit(file[grep('Tmax:', file)], '^Tmax: *')[[1]][2])
  
  # generation as factor
  # functions like max() can operate on ordered factors
  Data$gen <- factor(Data$gen, ordered=TRUE)
  
  # Number of traits (metagene handles either 1 or 2)
  # Try to guess the name of the file .001 where it reads the input parameters
  # otherwise, we try to infer it from the Data
#   fname <- strsplit(file[grep('NAME_OF_THIS_RUN', file)], c(' +'))[[1]][2]
  file1.path <- gsub('002', '001', file.path)
  if(file.exists(file1.path)){
    file1 = readLines(con=file1.path)
    n.traits <- as.integer(strsplit(file1[grep('Total_#_traits:', file1)],
                                    c(' +'))[[1]][2])
  } else {
    n.traits <- ifelse(with(Data, 
                            all(BV_Y==0) &
                            all(Dom_Y==0) &
                            all(phe_Y==0) &
                            all(indx_Y==0) &
                            all(indx_XY==0)),
                       1, 2)
  }
  if(n.traits==1) 
    Data <- Data[,-match(c('BV_Y', 'Dom_Y', 'phe_Y',
                                   'indx_Y', 'indx_XY'),
                                 names(Data))]
  
  # Number of individuals
  n.ind = max(Data$self)
  
  # Build object
  meta <- list(file.path     = file.path,
               file          = file,
               n.traits      = n.traits,
               n.generations = t.max,
               n.individuals = n.ind,
               Data          = Data)
  attr(meta, 'data.starting.line') <- data_start.idx
  class(meta) <- 'metagene'  
  return(meta)
}


#' @method summary metagene
#' @export
summary.metagene <- function(object, ...) {
#   attach(x)
  xsummary <- function(x) c(summary(x),
                            'SD'=sd(x),
                            'Var'=var(x))[c(4, 7, 8, 1:3, 5:6)]
  
  breeding.values <- list(
    global = do.call(rbind,
                     tapply(object$BV_X,
                            rep(1, object$n.individuals),
                            xsummary)),
    generation = do.call(rbind,
                         tapply(object$BV_X,
                                object$gen,
                                xsummary)),
    sex = do.call(rbind,
                  tapply(object$BV_X,
                                object$sex,
                                xsummary))
  )

  # Exclude founders except for generation-specific results
  object.noF <- object[object$gen!=0,]
  phenotype <- list(
    'global (exc. founders)' = do.call(rbind,
                                       tapply(object.noF$phe_X,
                                              rep(1, object.noF$n.individuals),
                                              xsummary)),
    generation = do.call(rbind,
                         tapply(object$phe_X,
                                object$gen,
                                xsummary)),
    'sex (exc. founders)' = do.call(rbind,
                                    tapply(object.noF$phe_X,
                                           object.noF$sex,
                                           xsummary))
  )
  
  # If only 1 generation, don't split by generation
  if(object$n.generations < 2) {
    breeding.values <-
      breeding.values[-which(names(breeding.values)=='generation')]
    phenotype <-
      phenotype[-which(names(phenotype)=='generation')]
  }
  
  # Has spatial structure?
  is.spatial <- exists('sp_X', where = as.environment(object$Data))
  
  summeta <- list(file.path       = object$file.path,
                  file            = object$file,
                  n.traits        = object$n.traits,
                  n.generations   = object$n.generations,
                  n.individuals   = object$n.individuals,
                  Data            = object$Data,
                  breeding.values = breeding.values,
                  phenotype       = phenotype,
                  is.spatial      = is.spatial)
  class(summeta) <- 'summary.metagene'  
  return(summeta)
}

#' @method print summary.metagene
#' @export
print.summary.metagene <- function(x, ...) {
  cat('Metagene simulated dataset\n===========================\n')
  cat('Number of individuals: ', format(x$n.individuals, width=4),  '\n')
  cat('Number of traits     : ', format(x$n.traits, width=4),       '\n')
  cat('Number of generations: ', format(x$n.generations, width=4),  '\n')
  cat('Spatial structure    : ', ifelse(x$is.spatial, 'Yes', 'No'), '\n')
  cat('Selection strategy: (Warning: fixed info)\t
      diagonal; 10+10 descendants per mating; select best 80+80\n')
  cat('\nBreeding values: ##########################\n')
  for(i in 1:length(x$breeding.values)){
    cat(names(x$breeding.values)[i], ':\n')
    print(x$breeding.values[[i]])
  }
  cat('\nPhenotypic values: ########################\n')
  for(i in 1:length(x$phenotype)){
    cat(names(x$phenotype)[i], ':\n')
    print(x$phenotype[[i]])
  }
  cat('\nHeritability: #############################\n')
#   browser()
  for(i in 1:length(x$phenotype)){
    cat(names(x$phenotype)[i], ':\n')
    if(is.null(dim(x$phenotype[[i]])))
      print(x$breeding.values[[i]][,'Var']/x$phenotype[[i]][,'Var'])
    else
      print(x$breeding.values[[i]][,'Var']/x$phenotype[[i]][,'Var'])
  }
}

#' Plot method for metagene objects
#' 
#' Plots either genetic and phenotypic values, or the spatial component of the
#' phenotype
#' @param x a metagene object.
#' @param type character. If 'default', the empirical density of the breeding and phenotypical values will be represented by generation. If 'spatial', the map of the spatial component will be plotted.
#' @param ... Further layers passed to \code{\link[ggplot2]{ggplot}}.
#' @method plot metagene
#' @import ggplot2
#' @export
plot.metagene <- function(x, type = c('default', 'spatial'), ...) {
#   dat <- data(x)
  type <- match.arg(type)
  
  if(type == 'spatial') {
    stopifnot('spatial' %in% names(x))
    spdat <- with(as.data.frame(x),
                  data.frame(irow, icol, z = sp_X,
                             model = 'spatial'))
    p <- spatial.plot(spdat, scale = 'div') + 
      facet_wrap(~ model)
  }
  else {
    dat <- data.frame(label=rep(c('genotype', 'phenotype'),
                                each = nindividuals(x)),
                      generation = rep(factor(x$gen), 2),
                      sex = rep(x$sex, 2),
                      value = c(x$BV_X, x$phe_X))
    p <- ggplot(dat, aes(x = value, fill = label)) +
      geom_density(alpha=.3) +
      facet_grid(generation~.) +
      labs(x = "Value by generation") 
  }
  
  if( !missing(...) ) {
    p <- p + ...
  }
  
  p
}

#### Interface functions ####

#' Extract the number of traits
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @export
get_ntraits <- function(x, ...) UseMethod('get_ntraits')
#' @export
get_ntraits.metagene <- function(x, ...) {
  return(x$n.traits)
}

#' Number of generations
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @export
ngenerations <- function(x, ...) UseMethod('ngenerations')
#' @export
ngenerations.metagene <- function(x, ...) {
  return(x$n.generations)
}

#' Number of individuals
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @export
nindividuals <- function(x, ...) UseMethod('nindividuals')
#' @export
nindividuals.metagene <- function(x, exclude.founders = FALSE, ...) {
  N <- x$n.individuals
  if(exclude.founders) N = N - sum(x$gen==0)
  return(N)
}

#' Get the Pedigree from an object
#' 
#' Returns an object from the formal class \code{pedigree}.
#' @param x object to extract pedigree from
#' @param ... Arguments to be passed to methods.
#' @references \code{\link[pedigreemm]{pedigree-class}} from package
#'   \code{pedigreemm}
#' @export
get_pedigree <- function(x, ...) UseMethod('get_pedigree')
#' @describeIn get_pedigree Get the pedigree from a metagene object
#' @export
get_pedigree.metagene <- function(x, ...) {
  return(with(x$Data, pedigreemm::pedigree(sire=dad, dam=mum, label=self)))
}

#' Coerce to a data.frame
#' 
#' This function returns a \code{\link[pedigree]{pedigree}} object in a data
#' frame
#' 
#' @method as.data.frame pedigree
#' @param x a \code{\link[pedigree]{pedigree}} object
#' @param ... not used
#'   
#' @return returns a data frame with one row per individual, the first column
#'   being the identification code, and the other two columns are dad and mum
#'   codes respectively.
#' @export
as.data.frame.pedigree <- function(x, ...) {
  y <- as(x, 'data.frame')
  codes <- as.numeric(row.names(y))
  z <- cbind(self=codes, sapply(y, function(x) codes[x]))
  return(as.data.frame(z))
}

#' Coerce to a data.frame
#' 
#' This function returns the data frame with the pedigree and phenotypes of the
#' simulated individuals.
#' 
#' @method as.data.frame metagene
#' @param x a metagene object
#' @param ... not used
#' @param exclude.founders logical: should the data.frame contain the genetic
#'   values of the founders?
#'   
#' @return returns a data frame with one row per individual, with the spatial
#'   coordinates if applicable, the pedigree information, the generation, the
#'   true breeding value, the phenotype, the sex, the spatially structured
#'   component of the phenotype and other internal metagene variables.
#' @export
as.data.frame.metagene <- function(x, ..., exclude.founders = TRUE) {
  # Exclude founders if appropriate
  if(exclude.founders) y <- x[x$gen!=0, ]
  else y <- x
  
  # If it has a spatial structure, add coordinates as columns
  if('spatial' %in% names(y))
    dat = data.frame(rbind(matrix(NA, sum(y$gen==0), 2), 
                           coordinates(y)), 
                     y$Data)
  else dat = y$Data
  
  return(dat)
}


#' Extract or replace data in a metagene object
#' 
#' @name Extract.metagene
#' @param x a metagene object.
#' @param name character. A varaible name.
#' @param value a vector.
#' @param ... a vector of integer indices or names of columns in the dataset.
NULL

#' Extract columns directly from the dataframe
#' @rdname Extract.metagene
#' @export
"$.metagene" <- function(x, name) {
  # Decide wether we asked for an element of the list
  # or an element of the dataframe
  if(name %in% names(x)) x[[name]]
  else x[['Data']][[name]]
}

#' Write columns of the dataframe
#' @rdname Extract.metagene
#' @export
"$<-.metagene" <- function(x, name, value) {
  # Decide wether we asked for an element of the list
  # or an element of the dataframe
  if(name %in% names(x)) x[[name]] <- value
  else x[['Data']][[name]] <- value
  x
}

#' Subset data
#' @rdname Extract.metagene
#' @export
"[.metagene" <- function(x, ...) {
  founders.idx <- which(x$gen==0)
  nothing = matrix(NA, max(founders.idx), 2)
  coord <- rbind(nothing, coordinates(x))
  Data.subset <- "[.data.frame"(cbind(coord, x$Data), ...)
  # update items
  y <- x
  y$Data <- Data.subset[, -(1:2)]
  y$n.generations <- max(y$gen)
  y$n.individuals <- nrow(Data.subset)
  # drop unused levels of factors
  y.factors <- which(sapply(y$Data, is.factor))
  for(var in y.factors) {
    y$Data[[var]] <- factor(y$Data[[var]], 
                            ordered=is.ordered(y$Data[[var]]))
  }
  # Keep coordinates only for the subset
  # (coordinates exclude founders)
  coord <- Data.subset[!apply(is.na(Data.subset[, 1:2]), 1, all), 1:2]
  if(nrow(coord) > 0)
    y$spatial$spatial.points <- sp::SpatialPoints(coord)
  else y$spatial$spatial.points <- list(NULL)
  return(y)
}

#' Breeding values
#' @param x a metagene object.
#' @export
b.values <- function(x) {
  stopifnot(inherits(x, 'metagene'))
  dat <- x$Data[,c('self', 'gen', 'sex', 'BV_X')]
  return(dat)
}

# # Fixed and random effects
# # Metagene simulates genetic selection, which  induces a generation effect
# # The variable sex is there, but there is no sex effect.
# # Apart from that, you can model the effects however you like.
# effects <- function(x) UseMethod('effects')
# effects.metagene <- function(x) {
#   return(list(generation = list(type = 'fixed', symbol = 'beta_g'),
#               sex = list(type = 'fixed', symbol = 'beta_s'),
#               genetic = list(type = 'random', symbol = 'u'))
#   )
# }


#### Simulate spatial structure ####

#' Simulate a spatial structure
#' 
#' Takes a metagene simulated dataset and distributes its individuals randomly 
#' into a more or less square spatial region. Furthermore, it takes part of the 
#' phenotypic noise and puts some spatial structure into it.
#' 
#' Founders are not put into place, as they don't have phenotypic values. 
#' Therefore, they are returned without coordinates nor spatial values.
#' 
#' The variance of the spatial field is given as a proportion of the variance of
#' the random noise that was added to the Breeding Values to produce the 
#' phenotypes. The phenotypes are afterwards recalculated to preserve the 
#' heritability.
#' 
#' The spatial unit is the distance between consecutive trees.
#' @param meta A metagene object
#' @param variance A number between 0 and 1. The variance of the spatial field 
#'   as a proportion of the non-inheritable pehnotype variance. See Details.
#' @param range A number between 0 and 1. The range of the spatial field, as a 
#'   proportion of the region size.
#' @param ... Arguments to be passed to methods.
#' @return Another metagene object with spatial structure given through an 
#'   additional column \code{sp_X} with the spatially structured component of 
#'   the phenotype, and a 'spatial' list element with coordinates and simulation
#'   information
#' @export
sim.spatial <- function(meta, variance, range, ...) UseMethod('sim.spatial')
#' @export
sim.spatial.metagene <- function(meta, variance = 0.5, range = 0.5, ...) {
  
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Pakage INLA needed for this simulating the spatial structure. Please install.",
         call. = FALSE)
  }
  
  #### Distribute the individuals (except founders) randomly in space
#   meta <- meta[meta$gen!=0,]
  founders.idx <- which(meta$gen==0)
  N <- nindividuals(meta) - length(founders.idx)
  # Try to split N into nrxnc, as square as possible
  nr <- floor(sqrt(N)); nc <- ceiling(N/nr)
  
  # Randomize individuals
  spatial.order <- sample(meta$self[-founders.idx], N)
  spatial.coord <- head(as.data.frame(INLA::inla.node2lattice.mapping(nr, nc)), N)
  SP <- sp::SpatialPoints(spatial.coord[order(spatial.order), ])
  
  #### Simulate spatial field
  
  # The spatial unit will be the distance between consecutive trees
  
  # Generate INLA mesh
  mesh <- INLA::inla.mesh.2d(loc = coordinates(SP),
                       cutoff = floor(.04*nr),
                       offset = floor(c(.05, .2)*nr),
                       max.edge = floor(c(.05, .1)*nr))
  # plot(mesh)
  
  # SPDE structure The variance of the spatial effect is a proportion of the 
  # phenotype's non-genetical variance (1-h) 
  BV.var <- var(meta$BV_X[meta$gen==1])
  ph.var <- var(meta$phe_X[meta$gen==1])
  h = BV.var/ph.var     # 1-h: proportion of phenotypical variance 
                        # due to non-genetical factors
  
  sigma0 = sqrt(variance * ph.var * (1-h))  ## field std.dev.
#   # Empirical correction: I know that the simulated field will have some
#   # border effect and that the variance in the area of interest will be lower.
#   sigma0 = 1.2 + sigma0
  range0 = range*nr    ## field range (a 'range' proportion of the region size)
  # Field parameters
  kappa0 = sqrt(8)/range0
  tau0 = 1/(sqrt(4*pi)*kappa0*sigma0)
  
  # SPDE object
  spde=INLA::inla.spde2.matern(mesh,
                         B.tau=cbind(log(tau0),1,0),
                         B.kappa=cbind(log(kappa0),0,1),
                         theta.prior.mean=c(0,0),
                         theta.prior.prec=1)
  
  # Precision matrix (using the prior means of the parameters)
  Q=suppressMessages(INLA::inla.spde2.precision(spde,theta=c(0,0)))
  
  # Sample
  x=suppressWarnings(as.vector(INLA::inla.qsample(n=1,Q)))
  # str(x)
  # The mesh vertices don't coincide with the observation locations
  # The A matrix provides a mapping such that y = Ax
#   A = INLA::inla.spde.make.A(mesh, loc = coordinates(SP))
  
#   # Visualize
#   plot(mesh, rgl=TRUE, col=x,
#        color.palette=function(n) grey.colors(n, 1, 0),
#        draw.edges=FALSE, draw.segments=TRUE, draw.vertices=FALSE)
#   
#   # Project into a 100x100 matrix (only inner area)
#   proj.mat = INLA::inla.mesh.projector(mesh, xlim=c(1, nr), ylim=c(1,nc))    
    # This contains the A matrix in proj$proj$A
#   ss.matrix <- INLA::inla.mesh.project(proj.mat, field = x)
#   image(ss.matrix)
  
  # Project into the observations locations
  proj.vec = INLA::inla.mesh.projector(mesh, loc=coordinates(SP))    
  ss.vec <- INLA::inla.mesh.project(proj.vec, field = x)
#   str(ss.vec)       # Value of the simulated spatial field in the obs locs
#   summary(ss.vec)
#   var(ss.vec)      # This should be similar  to BV.var
  
  #### Reduce variability from the phenotype and insert the spatial effect
  # to preserve the heritability
  spatial.var <- var(ss.vec)
  meta$sp_X <- c(NA[founders.idx], ss.vec)
  meta$phe_X[-founders.idx] <- 
    (meta$BV_X[-founders.idx] + ss.vec) + 
    sqrt((ph.var-BV.var-spatial.var)/(ph.var-BV.var))*
    (meta$phe_X - meta$BV_X)[-founders.idx]
  
  # # These two must be similar. i.e. kept total variance.
  # var(meta$phe_X)
  # var(meta$phe_sp_X)
  
  #### Return metagene object with spatial effect ####
  
  meta[['spatial']] <- list(spatial.points=SP, mesh=mesh, sim.field=x)
  return(meta)
}

#### Spatial interface functions ####

# sp::coordinates() is an S4 function
# Register the S3 class 'metagene' as an S4 class
setOldClass('metagene')
# @importFrom sp coordinates
# @importClassesFrom sp Spatial
# @importMethodsFrom sp coordinates coordinates<-
# @export
#' @import sp
setMethod('coordinates', signature = 'metagene', 
          function(obj, ...) {
            if(!('spatial' %in% names(obj)))
              stop("This metagene object has no spatial structure. 
                   Use sim.spatial().")
            coordinates(obj$spatial$spatial.points, ...)
          }
)

# Dummy method (Provide a proper def)
# @export
setMethod('coordinates<-', signature = 'metagene', 
          function(object, value) {
            NULL
          }
)
