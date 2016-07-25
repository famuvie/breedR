if (is.na(getOption("repos")["breedR"])) {
  expr <- expression(
    r <- getOption("repos"),
    if (is.null(r)) r["CRAN"] <- "https://cloud.r-project.org/",
    r["breedR"] <- "https://famuvie.github.io/breedR",
    options(repos = r)
  )
  
  ## set up repo for current and future sessions
  local(eval(expr))
  
  ## find out user's Rprofile
  rprofile <- ifelse(Sys.getenv("R_PROFILE_USER") == "",
                     file.path(Sys.getenv('HOME'), '.Rprofile'),
                     Sys.getenv("R_PROFILE_USER"))
  
  ## write expression
  cat(c("local({", sapply(expr, deparse), "})"),
        file = rprofile, sep = "\n", append = TRUE)
}
