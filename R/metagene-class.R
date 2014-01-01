#%%%%%%%%%%%%%%%%%%%%#
#%% metagene class %%#
#%% Facundo Muñoz  %%#
#%%%%%%%%%%%%%%%%%%%%#

# require(pedigreemm)   # For returning the pedigree in the pedigree::pedigree-class
# require(INLA)         # For simulating a spatial effect

#' Metagene Data Input
#' 
#' Read the output file (the one with the pedigree. Usually: insim.002)
#' of te Metagene program, and build an object of class metagene.
#' @param fname file name of the second metagene output file (usually: insim.002)
#' @return read.metagene returns the data in the file as an object of class \code{metagene}
#' @references \url{http://www.igv.fi.cnr.it/noveltree/}
#' @family metagene
#' @export
read.metagene <- function(fname) {
  file.path = fname
  file      = readLines(con=file.path)
  pedigree.start <- grep('pedigree', file) + 1
  pedigree <- read.table(file=file.path, quote="", skip=pedigree.start, col.names=scan(file=file.path, what=character(), skip=pedigree.start - 1, nlines=1, sep=','))
  t.max = as.integer(strsplit(file[grep('Tmax:', file)], '^Tmax: *')[[1]][2])
  
  # generation as factor
  # functions like max() can operate on ordered factors
  pedigree$gen <- factor(pedigree$gen, ordered=TRUE)
  
  # Number of traits (metagene handles either 1 or 2)
  # We try to guess the name of the file .001 where it reads the input parameters
  # otherwise, we try to infer it from the pedigree
#   fname <- strsplit(file[grep('NAME_OF_THIS_RUN', file)], c(' +'))[[1]][2]
  file1.path <- gsub('002', '001', file.path)
  if(file.exists(file1.path)){
    file1 = readLines(con=file1.path)
    n.traits <- as.integer(strsplit(file1[grep('Total_#_traits:', file1)], c(' +'))[[1]][2])
  } else {
    n.traits <- ifelse(with(pedigree, all(BV_Y==0) & all(Dom_Y==0) & all(phe_Y==0) & all(indx_Y==0) & all(indx_XY==0)), 1, 2)
  }
  if(n.traits==1) 
    pedigree <- pedigree[,-match(c('BV_Y', 'Dom_Y', 'phe_Y',
                                   'indx_Y', 'indx_XY'),
                                 names(pedigree))]
  
  # Number of individuals
  n.ind = max(pedigree$self)
  
  # Build object
  meta <- list(file.path=file.path, file=file, n.traits=n.traits, n.generations=t.max, n.individuals=n.ind, pedigree=pedigree)
  attr(meta, 'pedigree.starting.line') <- pedigree.start
  class(meta) <- 'metagene'  
  return(meta)
}


#' Summary method for metagene objects
#' 
#' Prints a summary of a metagene objects.
#' @method summary metagene
#' @param x A metagene object
#' @return Prints summary
#' @export
summary.metagene <- function(object, ...) {
#   attach(x)
  xsummary <- function(x) c(summary(x), 'SD'=sd(x), 'Var'=var(x))[c(4, 7, 8, 1:3, 5:6)]
  
  breeding.values <- list(
    global = do.call(rbind, tapply(object$BV_X, rep(1, object$n.individuals), xsummary)),
    generation = do.call(rbind, tapply(object$BV_X, object$gen, xsummary)),
    sex = do.call(rbind, tapply(object$BV_X, object$sex, xsummary))
  )

  # Exclude founders except for generation-specific results
  object.noF <- object[object$gen!=0,]
  phenotype <- list(
    'global (exc. founders)' = do.call(rbind, tapply(object.noF$phe_X, rep(1, object.noF$n.individuals), xsummary)),
    generation = do.call(rbind, tapply(object$phe_X, object$gen, xsummary)),
    'sex (exc. founders)' = do.call(rbind, tapply(object.noF$phe_X, object.noF$sex, xsummary))
  )
  
  # If only 1 generation, don't split by generation
  if(object$n.generations < 2) {
    breeding.values <- breeding.values[-which(names(breeding.values)=='generation')]
    phenotype <- phenotype[-which(names(phenotype)=='generation')]
  }
  
  summeta <- list(file.path=object$file.path, file=object$file, n.traits=object$n.traits, n.generations=object$n.generations, n.individuals=object$n.individuals, pedigree=object$pedigree, breeding.values=breeding.values, phenotype=phenotype)
  class(summeta) <- 'summary.metagene'  
  return(summeta)
}

#' @export
print.summary.metagene <- function(x, ...) {
  cat('Metagene simulated dataset\n===========================\n')
  cat('Number of individuals: ', format(x$n.individuals, width=4), '\n')
  cat('Number of traits     : ', format(x$n.traits, width=4), '\n')
  cat('Number of generations: ', format(x$n.generations, width=4), '\n')
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
#' Plots either genetic and phenotypic values, or the spatial component of the phenotype
#' @method plot metagene
#' @export
plot.metagene <- function(x, type = c('default', 'spatial'), ...) {
#   dat <- data(x)
  type <- match.arg(type)
  n = nindividuals(x)
  
  if(type == 'spatial') {
    stopifnot('spatial' %in% names(x))
    ggplot(as.data.frame(x), aes(icol, irow)) + 
      geom_tile(aes(fill=sp_X)) + 
      scale_fill_gradient(low='green', high='red')
  }
  else {
    dat <- data.frame(label=rep(c('genotype', 'phenotype'), each=n),
                      generation = rep(factor(x$gen), 2),
                      sex = rep(x$sex, 2),
                      value = c(x$BV_X, x$phe_X))
    ggplot(dat, aes(x = value, fill = label)) + geom_density(alpha=.3) + facet_grid(generation~.) + labs(x = "Value by generation") 
  }
}

#### Interface functions ####

#' Number of traits
#' @export
ntraits <- function(x) UseMethod('ntraits')
#' @export
ntraits.metagene <- function(x) {
  return(x$n.traits)
}

#' Number of generations
#' @export
ngenerations <- function(x) UseMethod('ngenerations')
#' @export
ngenerations.metagene <- function(x) {
  return(x$n.generations)
}

#' Number of individuals
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
#' Returns an object from the formal class 'pedigree'
#' @export
get_pedigree <- function(x, ...) UseMethod('get_pedigree')
get_pedigree.metagene <- function(x, ...) {
  return(with(x$pedigree, pedigreemm::pedigree(sire=dad, dam=mum, label=self)))
}

#' Coerce to a data.frame
#' 
#' This function returns a \code{\link[pedigree]{pedigree}} object in a data frame
#' 
#' @method as.data.frame metagene
#' @param x a \code{\link[pedigree]{pedigree}} object
#' @param ... not used
#' 
#' @return returns a data frame with one row per individual, the first column being the identification code, and the other two columns are dad and mum codes respectively.
#' @export
as.data.frame.pedigree <- function(x, ...) {
  y <- as(x, 'data.frame')
  z <- cbind(self=as.numeric(row.names(y)), y)
}

# # data (method masked from utils)
# # Returns a dataframe
# # this causes trouble. Better use the "[" method.
# data <- function(x, ...) UseMethod('data')
# data.metagene <- function(x) {
#   return(x$pedigree)
# }
# data.default <- function(...) {
#   utils::data(...)
# }

#' Coerce to a data.frame
#' 
#' This function returns the data frame with the pedigree and phenotypes of the simulated individuals.
#' 
#' @method as.data.frame metagene
#' @param x a metagene object
#' @param ... not used
#' @param exclude.founders logical: should the data.frame contain the genetic values of the founders?
#' 
#' @return returns a data frame with one row per individual, with the spatial coordinates if applicable, the pedigree information, the generation, the true breeding value, the phenotype, the sex, the spatially structured component of the phenotype and other internal metagene variables.
#' @export
as.data.frame.metagene <- function(x, ..., exclude.founders = TRUE) {
  # Exclude founders if appropriate
  if(exclude.founders) x = x[x$gen!=0, ]
  
  # If it has a spatial structure, add coordinates as columns
  if('spatial' %in% names(x))
    dat = data.frame(rbind(matrix(NA, sum(x$gen==0), 2), 
                           coordinates(x)), 
                     x$pedigree)
  else dat = x$pedigree
  
  return(dat)
}


#' Extract columns directly from the dataframe
#' 
#' @rdname Extract.metagene
#' @export
"$.metagene" <- function(x, name) {
  # Decide wether we asked for an element of the list
  # or an element of the dataframe
  if(name %in% names(x)) x[[name]]
  else x[['pedigree']][[name]]
}

#' Write columns of the dataframe
#' 
#' @rdname Extract.metagene
#' @export
"$<-.metagene" <- function(x, name, value) {
  # Decide wether we asked for an element of the list
  # or an element of the dataframe
  if(name %in% names(x)) x[[name]] <- value
  else x[['pedigree']][[name]] <- value
  x
}

#' Subset data
#' 
#' @rdname Extract.metagene
#' @export
"[.metagene" <- function(x, ...) {
  pedigree.subset <- "[.data.frame"(x$pedigree, ...)
  # update items
  y <- x
  y$pedigree <- pedigree.subset
  y$n.generations <- max(y$gen)
  y$n.individuals <- nrow(pedigree.subset)
  # drop unused levels of factors
  y.factors <- which(sapply(pedigree.subset, is.factor))
  for(var in y.factors) {
    y$pedigree[[var]] <- factor(pedigree.subset[[var]], ordered=is.ordered(pedigree.subset[[var]]))
  }
  return(y)
}

#' Breeding values
#' @export
b.values <- function(x, ...) {
  stopifnot(inherits(x, 'metagene'))
  dat <- data(x)[,c('self', 'gen', 'sex', 'BV_X')]
}

# # Fixed and random effects
# # Metagene simulates genetic effects, which in turn induces a generation effect
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
sim.spatial <- function(meta, ...) UseMethod('sim.spatial')
#' @method sim.spatial metagene
#' @export
sim.spatial.metagene <- function(meta, ...) {
  #### Distribute the individuals (except founders) randomly in space
#   meta <- meta[meta$gen!=0,]
  founders.idx <- which(meta$gen==0)
  N <- nindividuals(meta) - length(founders.idx)
  # Try to split N into nrxnc, as square as possible
  nr <- floor(sqrt(N)); nc <- ceiling(N/nr)
  
  # Randomize individuals
  spatial.order <- sample(meta$self[-founders.idx], N)
  spatial.coord <- head(as.data.frame(inla.node2lattice.mapping(nr, nc)), N)
  SP <- SpatialPoints(spatial.coord[order(spatial.order), ])
  
  #### Simulate spatial field
  
  # The spatial unit will be the distance between consecutive trees
  
  # Generate INLA mesh
  mesh <- inla.mesh.2d(loc = coordinates(SP),
                       cutoff = floor(.04*nr),
                       offset = floor(c(.05, .2)*nr),
                       max.edge = floor(c(.05, .1)*nr))
  # plot(mesh)
  
  # SPDE structure
  # We give the spatial effect a variance of 30% of phenotype's variance 
  # (the same as the heritability. 
  # i.e. proportion between the genetic variance and phenotype's variance)
  # But we reduce the phenotypic variance in the same amount, 
  # to preserve the heritability
  BV.var <- var(meta$BV_X[meta$gen==1])
  ph.var <- var(meta$phe_X[meta$gen==1])
  h = BV.var/ph.var     # 1-h: proportion of phenotypical variance due to environmental factors
  
  sigma0 = sqrt(BV.var)  ## field std.dev. (of the same order of magnitude than genetic variance)
  range0 = 0.2*nr        ## field range (about 20% of the region size)
  # Field parameters
  kappa0 = sqrt(8)/range0
  tau0 = 1/(sqrt(4*pi)*kappa0*sigma0)
  
  # SPDE object
  spde=inla.spde2.matern(mesh,
                         B.tau=cbind(log(tau0),1,0),
                         B.kappa=cbind(log(kappa0),0,1),
                         theta.prior.mean=c(0,0),
                         theta.prior.prec=1)
  
  # Precision matrix (using the prior means of the parameters)
  Q=suppressMessages(inla.spde2.precision(spde,theta=c(0,0)))
  
  # Sample
  x=suppressWarnings(as.vector(inla.qsample(n=1,Q)))
  # str(x)
  # The mesh vertices don't coincide with the observation locations
  # The A matrix provides a mapping such that y = Ax
#   A = inla.spde.make.A(mesh, loc = coordinates(SP))
  
#   # Visualize
#   plot(mesh, rgl=TRUE, col=x,
#        color.palette=function(n) grey.colors(n, 1, 0),
#        draw.edges=FALSE, draw.segments=TRUE, draw.vertices=FALSE)
#   
#   # Project into a 100x100 matrix (only inner area)
#   proj.mat = inla.mesh.projector(mesh, xlim=c(1, nr), ylim=c(1,nc))    # This contains the A matrix in proj$proj$A
#   ss.matrix <- inla.mesh.project(proj.mat, field = x)
#   image(ss.matrix)
  
  # Project into the observations locations
  proj.vec = inla.mesh.projector(mesh, loc=coordinates(SP))    
  ss.vec <- inla.mesh.project(proj.vec, field = x)
#   str(ss.vec)       # Value of the simulated spatial field in the observation locations
#   summary(ss.vec)
#   var(ss.vec)      # This should be similar  to BV.var
  
  #### Reduce variability from the phenotype and insert the spatial effect
  spatial.var <- var(ss.vec)
  meta$sp_X <- c(NA[founders.idx], ss.vec)
  meta$phe_X[-founders.idx] <- (meta$BV_X[-founders.idx] + ss.vec) + sqrt((ph.var-BV.var-spatial.var)/(ph.var-BV.var))*(meta$phe_X - meta$BV_X)[-founders.idx]
  
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
setMethod('coordinates', signature = 'metagene', 
          function(obj, ...) {
            if(!('spatial' %in% names(obj)))
              stop("This metagene object has no spatial structure. Use sim.spatial().")
            coordinates(obj$spatial$spatial.points, ...)
          }
)

