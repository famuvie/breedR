old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

#### Context: Models with several effects working together ####
context("Models with several effects working together") 

# TODO:
# mixed effects models * spatial models * genetic effects