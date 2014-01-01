### Functions intended to aid in the development process
### and in the documentation of the package

#' Plots a flowchart of function calls
#' 
#' This function is a wrapper for mvbutils::foodweb()
funcall <- function(fun) {
  mvbutils::foodweb(where='package:breedR', prune = fun,
                    border = TRUE, expand.xbox = 1.2,
                    boxcolor = "#FC6512", textcolor = "black",
                    cex = 1, lwd=2)
}