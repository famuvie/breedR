---
layout: default
title: Home
---

This package provides frequentist and Bayesian statistical tools to build predictive models useful for the breeders, quantitative genetists and forest genetic resources analysts communities. It aims to assess the genetic value of individuals under a number of situations, including spatial autocorrelation, genetic/environment interaction and competition. It is under active development as part of the [Trees4Future project](http://www.trees4future.eu/ "T4F"), particularly developed having forest genetic trials in mind. But can be used for animals or other situations as well.

Next [Training Workshop](workshop)

### Installation
This will install the latest stable version of breedR.
For the latest development version, please refer to the [development site](https://github.com/famuvie/breedR)

**Requirements:**
- Install and get started with [R](getR)
- Install the `devtools` R-package: `install.packages('devtools')`

{% highlight r %}
## Installation of breedR (latest stable version)
library(devtools)
install_github('famuvie/breedR', ref = github_release())
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

#### Citing
- If you use this package please cite it
- `citation('breedR')`
