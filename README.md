breedR
======

### Statistical methods for forest genetic resources analysts

This package provides frequentist and Bayesian statistical tools to build predictive models useful for the breeders, quantitative genetists and forest genetic resources analysts communities. It aims to assess the genetic value of individuals under a number of situations, including spatial autocorrelation, genetic/environment interaction and competition. It is under active development as part of the [Trees4Future project](http://www.trees4future.eu/ "T4F"), particularly developed having forest genetic trials in mind. But can be used for animals or other situations as well.

If you have questions, please join our discussion [group](http://groups.google.com/group/breedr)

#### Usage and collaboration

##### First time only
- Install [R](http://cran.r-project.org/ "CRAN") from CRAN
- If you work under MS Windows platform
  - Install [`RTools`](http://cran.r-project.org/bin/windows/Rtools/) as well
- Install the `devtools` R-package
  - `install.packages('devtools')`
- For beta-testers;
  - Create a [GitHub](https://github.com/join) account
  - Install the `testthat` R-package
    - `install.pacakges('testthat')`

##### Installation of breedR (latest dev. version)
  - `devtools::install_github('famuvie/breedR')`

##### How to use
  ```R
  library(breedR)             # Load the package
  news(package = 'breedR')    # Check the changelog
  example('breedR')           # Check-out the basic example
  demo(package='breedR')      # Available demos on features
  demo('Metagene-spatial')    # Execute some demos
  demo('globulus')
  ```

##### Test cycle
- Check the automated tests
  ```R
  library('testthat')
  test_package('breedR')
  ```
  
- Try it with your own data or with provided datasets
- Report [issues](https://github.com/famuvie/breedR/issues "Issues page")

##### Citing
- If you use this package please cite it
- `citation('breedR')`
