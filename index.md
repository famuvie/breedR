---
layout: default
title: breedR
---

This package provides frequentist and Bayesian statistical tools to build predictive models useful for the breeders, quantitative genetists and forest genetic resources analysts communities. It aims to assess the genetic value of individuals under a number of situations, including spatial autocorrelation, genetic/environment interaction and competition. It is under active development as part of the [Trees4Future](http://www.trees4future.eu/ "T4F") and [ProCoGen](http://www.procogen.eu/) projects, particularly developed having forest genetic trials in mind. But can be used for animals or other situations as well.

### Installation
This will install the latest stable version of breedR. For the latest development version, please refer to the [development site](https://github.com/famuvie/breedR)

-   Install and get started with [R](getR)
-   Set up breedR repository. This is only needed the first time.

{% highlight r %}
expr <- expression(
          r <- getOption("repos"),
          r["breedR"] <- "http://famuvie.github.io/breedR",
          options(repos = r)
        )
eval(expr)
rprofile <- file.path(Sys.getenv('HOME'), '.Rprofile')
if (file.exists(rprofile)) {
  p <- readLines(rprofile)
} else p <- c()
writeLines(c(p, sapply(expr, deparse)), rprofile)

{% endhighlight r %}

-   Install as a regular package. You can use package-management menus on RGui, Rstudio or whatever, or in plain R:
  
    {% highlight r %}
    install.packages('breedR')
    {% endhighlight r %}


### Getting started
- Check the [breedR-wiki](https://github.com/famuvie/breedR/wiki)
- ... particularly, the [Overview](https://github.com/famuvie/breedR/wiki/Overview) ([pdf](https://github.com/famuvie/breedR/wiki/Overview.pdf))

{% highlight r %}
library(breedR)             # Load the package
news(package = 'breedR')    # Check the changelog
example('breedR')           # Check-out the basic example
demo(package='breedR')      # Available demos on features
demo('Metagene-spatial')    # Execute some demos
demo('globulus')
{% endhighlight %}

- Any question? join our discussion [group](http://groups.google.com/group/breedr)
- Found a bug? help us by filing an [issue](https://github.com/famuvie/breedR/issues "Issues page")


### Training workshop materials

- [Poland Local Training Workshop](workshop_IBL) (December 2015, Warsaw) 
- [International Official Training Workshop](workshop) (June 2015, Jaca, Spain)


### Citing
- If you use this package please cite it
- `citation('breedR')`
