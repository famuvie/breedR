breedR
======

### Statistical methods for forest genetic resources analysts.

#### Description
This package provides frequentist and Bayesian statistical tools to build predictive models useful for the breeders, quantitative genetists and forest genetic resources analysts communities. It aims to assess the genetic value of individuals under a number of situations, including spatial autocorrelation, genetic/environment interaction and competition. It is under active development as part of the [Trees4Future project](http://www.trees4future.eu/ "T4F"), particularly developed having forest genetic trials in mind. But can be used for animals or other situations as well.

#### Usage and collaboration

##### First time only
- Install [R](http://cran.r-project.org/ "CRAN") from CRAN
- If you work under MS Windows platform
  - Install [`RTools`](http://cran.r-project.org/bin/windows/Rtools/) as well
- Install the `devtools` R-package
  - `install.packages('devtools')`
- Update to the latest development version
  - (Not estrictly necessary under Windows)
  - `devtools::install_github('devtools')`
- For beta-testers;
  - Create a [GitHub](https://github.com/join) account
  - Install the `testthat` R-package
    - `install.pacakges('testthat')

##### Installation
  - `library(devtools)`
  - `install_github('breedR', 'famuvie')`

##### How to use
- Load the package
  - `library(breedR)`
- Check the changelog
  - `news(package = 'breedR')`
- Check-out the basic example
  - `example('breedR')`
- Run some of the `demo(package='breedR')`
  - e.g. `demo('Metagene-spatial')`
  - or,  `demo('globulus')`

##### Test cycle
- Check the automated tests
  - `library('testthat')`
  - `test_dir(system.file('test', package = 'breedR'))`
- Try it with your own data or with provided datasets
- Report [issues](https://github.com/famuvie/breedR/issues "Issues page")

##### Citing
- If you use this package please cite it
- `citation('breedR')`