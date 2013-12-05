#' Metagene synthesized dataset
#' 
#' Simulated progeny of four generations under phenotypic selection and spatial arrangement.
#' The dataset includes the true Breeding values, the environmental effect and the observed phenotype for the 5600 individuals.
#' 
#' The generation of founders (generation 0) consist in 80 independent couples, each of which produces 20 descendants (10 of each sex).
#' The second generation descend from 80 random couples of the best individuals among the 1600 members of the first generation.
#' The same procedure is simulated to obtain the third and fourth generations.
#' 
#' The full dataset contains data for the \eqn{2\times 80 + 4\times 1600 = 6560}{2*80 + 4*1600 = 6560} simulated individuals.
#' However, the Metagene program does not simulate phenotypes for the founders.
#' 
#' The 6400 individuals from generations 1 to 4 were arranged at random in a \eqn{80\times 80}{80 x 80} spatial grid
#' 
#' 
#' 
#' @name m4
#' @docType data
#' @examples
#'   # Load, summarize and visualize data
#'   data(m4)
#'   summary(m4)
#'   plot(m4)
#'   
#'   # Environmental component of the phenotype
#'   # (spatially structured)
#'   plot(m4, type = 'spatial')
NULL