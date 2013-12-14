breedR
======

### Statistical methods for forest genetic resources analysts.
### Frequentist and Bayesian methods for breeders, quantitative genetists and forest genetic resources analysts.

#### Description
This package provides frequentist and Bayesian statistical tools to build predictive models useful for the breeders, quantitative genetists and forest genetic resources analysts communities. It aims to assess the genetic value of individuals under a number of situations, including spatial autocorrelation, genetic/environment interaction and competition. It is under active development as part of the Trees4Future project, particularly developed having forest genetic trials in mind. But can be used for animals or other situations as well.

#### Usage and collaboration

##### First time only
- Install [R](http://cran.r-project.org/ "CRAN")
  - (If you work under MS Windows platform, you need `RTools` in addition to the `base` distribution)
- Create a [GitHub](https://github.com/join) account
  - (If you don't have one, and want to provide feedback)
- Install the `devtools` R-package
  - `install.packages('devtools')`
- Update to the latest development version
  - (Not estrictly necessary under Windows)
  - `devtools::install_github('devtools')`

##### Test cycle
- Load `devtools`
  - `library(devtools)`
- Install or update the `breedR` package
  - `install_github('breedR', 'famuvie')`
- Check-out the basic example
  - `example('breedR', 'breedR')`
- Try it with your own data or with provided datasets
- Report [issues](https://github.com/famuvie/breedR/issues "Issues page")
