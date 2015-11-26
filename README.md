[![Travis-CI Build Status](https://travis-ci.org/famuvie/breedR.png?branch=master)](https://travis-ci.org/famuvie/breedR)
[![DOI](https://zenodo.org/badge/4357/famuvie/breedR.svg)](https://zenodo.org/badge/latestdoi/4357/famuvie/breedR)

breedR
======

### Statistical methods for forest genetic resources analysts

This [R](http://cran.r-project.org/ "CRAN") package provides frequentist and
Bayesian statistical tools to build predictive models useful for the breeders,
quantitative genetists and forest genetic resources analysts communities. It
aims to assess the genetic value of individuals under a number of situations,
including spatial autocorrelation, genetic/environment interaction and
competition. It is under active development as part of the [Trees4Future
project](http://www.trees4future.eu/ "T4F"), particularly developed having
forest genetic trials in mind. But can be used for animals or other situations
as well.

If you have questions, please join our discussion
[group](http://groups.google.com/group/breedr)

This site is concerned with the development and testing of breedR. If you want 
to use the most stable version of breedR, please check the [dissemination
site](http://famuvie.github.io/breedR/).

#### Installation

This will install the latest development version of `breedR`. Note that updating
requires explicitly repeating this operation. For regular use and management of
the package, follow the installation instructions at the [dissemination
site](http://famuvie.github.io/breedR/).

```R
devtools::install_github('famuvie/breedR')
```

#### Getting started
Check the [breedR-wiki](https://github.com/famuvie/breedR/wiki)
```R
library(breedR)                      # Load the package
browseVignettes(package = 'breedR')  # Read tutorials
news(package = 'breedR')             # Check the changelog
example('breedR')                    # Check-out the basic example
demo(package='breedR')               # Available demos on features
demo('Metagene-spatial')             # Execute some demos
demo('globulus')
```

#### Test cycle
breedR is in [beta](https://en.wikipedia.org/wiki/Development_stage#Beta) stage. Collaboration is welcome!
- Check the automated tests
    ```R
    library('testthat')
    test_package('breedR')
    ```
  
- Try it with your own data or with provided datasets
- Report [issues](https://github.com/famuvie/breedR/issues "Issues page")

#### Citing
- If you use this package please cite it
- `citation('breedR')`
