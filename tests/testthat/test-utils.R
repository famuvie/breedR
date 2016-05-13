### Test the auxiliar functions in utils.R ###

old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

context("Auxiliar functions")

