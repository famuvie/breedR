#' Longitudinal Larix dataset
#' 
#' Repeated measurements of microdensity data along 16 years, with climatic
#' covariates.
#' 
#' Each one of the 8 replicate is composed of 5 incomplete blocks, which gives a
#' nested variable with 40 levels.
#' 
#' @format A dataframe with 11897 observations of the followig variables
#' \itemize{
#'   \item{\code{self}}{   id of the tree}
#'   \item{\code{dad}}{    id of sire}
#'   \item{\code{mum}}{    id of dam}
#'   \item{\code{x, y}}{   coordinates (in rows/cols)}
#'   \item{\code{rep}}{    factor replicate with 8 levels}
#'   \item{\code{bl}}{     factor block with 40 levels nested within rep}
#'   \item{\code{yr}}{     factor (growth) year with 16 ordered levels}
#'   \item{\code{map}}{    mean annual precipitation}
#'   \item{\code{mat}}{    mean annual temperature}
#'   \item{\code{mi}}{     martone index}
#'   \item{\code{LAS}}{    phenotype LAS}
#'   \item{\code{DOS}}{    phenotype DOS}
#' }
#' 
#' @name larix
#' @docType data
#' @examples 
#' library(tidyr)
#' library(dplyr)
#' library(ggplot2)
#' data(larix)
#' 
#' ## N observations by year and replicate
#' with(larix, table(yr, rep))
#' 
#' ## Mean response evolution by replicate
#' larix %>%
#'   group_by(yr, rep) %>% 
#'   summarise(MLAS = mean(LAS)) %>% 
#'   ggplot(aes(yr, MLAS, group = rep)) +
#'   geom_line()
#'  
#' ## Visualise trial by year
#' ggplot(larix, aes(x, y)) +
#'   geom_tile(aes(fill = LAS)) + 
#'   facet_wrap(~ yr)
#' 
#' ## Correlations with environmental variables
#' if (require(GGally)) {
#'   ggpairs(larix[, c('map', 'mat', 'mi', 'LAS', 'DOS')])
#' }
NULL


#' Multi-site Douglas-fir dataset
#' 
#' Provenance test with measurements of height and circumference
#' 
#' The dataset includes measurements of height, circumference, angle and
#' branching, for a progeny from 11 different origins which were planted in 1998
#' in 3 different sites.
#' 
#' Only 4 out of 11 origins are available in all the three sites.
#' 
#' Circumference was measured in all the three sites in 2013, but 
#' measurements of height were taken in 2002, 2003 and 2004 for site 3, and in 
#' 2005 for site 2. Measurements of angle and  branching are only available for
#' site 3.
#' 
#' The spacing between trees is different among sites, and in general between
#' rows and columns. Coordinate units are meters.
#' 
#' @format A dataframe with 1021 individuals and the following 9 variables
#' \itemize{
#'   \item{\code{self}}{   id of the tree}
#'   \item{\code{dad}}{    id of sire or 0 if unknown}
#'   \item{\code{mum}}{    id of dam or 0 if unknown}
#'   \item{\code{orig}}{   factor origin of the seeds with 11 levels
#'    partially crossed with site}
#'   \item{\code{site}}{   factor site with 3 levels}
#'   \item{\code{block}}{  factor block with 127 levels nested within site}
#'   \item{\code{x, y}}{   coordinates (in m)}
#'   \item{\code{H02:H05}}{heights as measured in 2002-2005}
#'   \item{\code{C13}}{    circumference as measured in 2013}
#'   \item{\code{AN}}{     factor angle quality with 5 levels}
#'   \item{\code{BR}}{     factor branching quality with 5 levels}
#' }
#' 
#' @name douglas
#' @docType data
#' @examples 
#' data(douglas)
#' 
#' ## Individuals by origin and site
#' with(douglas, table(orig, site))
#' 
#' ## Measurements by variable and site
#' library(tidyr)
#' library(dplyr)
#' douglas %>%
#'  gather(variable, value, H02:BR, na.rm = TRUE) %>% 
#'  with(., table(variable, site))
#'  
#' ## Visualise trials
#' library(ggplot2)
#' ggplot(douglas, aes(x, y)) +
#'   geom_point() + 
#'   facet_wrap(~ site)
#'  
NULL

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
#'   library(ggplot2)
#'   qplot(BV_X, phe_X-BV_X, colour = dad, data = as.data.frame(m1)) + 
#'     geom_abline(intercept=0, slope=0, col='gray')
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
#'   library(ggplot2)
#'   qplot(BV_X, phe_X-BV_X, facets = .~gen, data = as.data.frame(m4)) + 
#'     geom_abline(intercept=0, slope=0, col='gray')
NULL
