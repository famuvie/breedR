### Test the computation of model matrices ###

old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

context("Model Matrix")

## TODO

## Test also the auxiliar function matrix.short16
## used to convert format from 16 cols to sparse.
