#' Eucalyptus Globulus dataset
#' 
#' Open-polinated field test with one generation of 1021 individuals.
#' 
#' The individuals are split into 14 genetic groups and arranged in blocks.
#' The plantation is gridded with a separation of 3 m.
#' 
#' @format A dataframe with 1021 individuals and the following 9 variables
#' \itemize{
#'   \item{\code{self}}{ id of the tree}
#'   \item{\code{dad}}{  id of sire or 0 if unknown}
#'   \item{\code{mum}}{  id of dam or 0 if unknown}
#'   \item{\code{gen}}{  generation (there is only 1)}
#'   \item{\code{gg}}{   genetic group}
#'   \item{\code{bl}}{   block}
#'   \item{\code{phe_X}}{observed phenotype}
#'   \item{\code{x, y}}{ coordinates (in m)}
#' }
#' @name globulus
#' @docType data
#' @examples data(globulus)
NULL

#' A small Metagene synthesized dataset
#' 
#' Simulated progeny of one single generation under phenotypic selection and spatial arrangement.
#' The dataset includes the true Breeding values, the environmental effect and the observed phenotype for 1760 individuals.
#' 
#' The generation of founders (generation 0) consist in 80 independent couples, each of which produces 20 descendants (10 of each sex).
#' 
#' The full dataset contains data for the \eqn{2\times 80 + 1600 = 1760}{2*80 + 1600 = 1760} simulated individuals.
#' However, the Metagene program does not simulate phenotypes for the founders.
#' 
#' The 1600 descendants were arranged at random in a \eqn{40\times 40}{40 x 40} spatial grid
#' 
#' 
#' 
#' @name m1
#' @docType data
#' @examples
#'   # Load, summarize and visualize data
#'   data(m1)
#'   summary(m1)
#'   plot(m1)
#'   
#'   # Environmental component of the phenotype
#'   # (spatially structured)
#'   plot(m1, type = 'spatial')
#'   
#'   # The phenotypes are noisy observations around the 
#'   # true Breeding values with a standard deviation of about 7.3.
#'   qplot(BV_X, phe_X-BV_X, colour = dad, data = as.data.frame(m1)) + 
#'     geom_abline(int=0, sl=0, col='gray')
NULL

#' Metagene synthesized dataset with four generations
#' 
#' Simulated progeny of four generations under phenotypic selection and spatial arrangement.
#' The dataset includes the true Breeding values, the environmental effect and the observed phenotype for the 6560 individuals.
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
#'   
#'   # The phenotypes are noisy observations around the Breeding values
#'   # with a standard deviation of about 10.
#'   qplot(BV_X, phe_X-BV_X, facets = .~gen, data = as.data.frame(m4)) + 
#'     geom_abline(int=0, sl=0, col='gray')
NULL