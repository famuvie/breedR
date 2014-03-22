\docType{data}
\name{m4}
\alias{m4}
\title{Metagene synthesized dataset with four generations}
\description{
  Simulated progeny of four generations under phenotypic
  selection and spatial arrangement. The dataset includes
  the true Breeding values, the environmental effect and
  the observed phenotype for the 6560 individuals.
}
\details{
  The generation of founders (generation 0) consist in 80
  independent couples, each of which produces 20
  descendants (10 of each sex). The second generation
  descend from 80 random couples of the best individuals
  among the 1600 members of the first generation. The same
  procedure is simulated to obtain the third and fourth
  generations.

  The full dataset contains data for the \eqn{2\times 80 +
  4\times 1600 = 6560}{2*80 + 4*1600 = 6560} simulated
  individuals. However, the Metagene program does not
  simulate phenotypes for the founders.

  The 6400 individuals from generations 1 to 4 were
  arranged at random in a \eqn{80\times 80}{80 x 80}
  spatial grid
}
\examples{
# Load, summarize and visualize data
  data(m4)
  summary(m4)
  plot(m4)

  # Environmental component of the phenotype
  # (spatially structured)
  plot(m4, type = 'spatial')

  # The phenotypes are noisy observations around the Breeding values
  # with a standard deviation of about 10.
  qplot(BV_X, phe_X-BV_X, facets = .~gen, data = as.data.frame(m4)) +
    geom_abline(int=0, sl=0, col='gray')
}
\keyword{datasets}

