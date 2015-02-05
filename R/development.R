## Functions for the development of breedR
## Not exported, and not even packaged 
## i.e. this file is listed in .Rbuildignore

## Render vignettes
## formats: pdf for the package
##          md  for the GitHub wiki
breedR.update_guides <- function(pkg = "."){
  pkg <- as.package(pkg)
  vigns <- tools::pkgVignettes(dir = pkg$path)

  ## Check existence of any vignette
  if (length(vigns$doc) == 0) 
    return()
  
  ## Otherwise, go ahead
  message("Building ", pkg$package, " vignettes")
  out <- breedR.build_vignettes(vigns)
  
  wiki_dir <- file.path(pkg$path, '..', 'breedR-wiki')
  breedR.move_vignettes(pkg, vigns, out, wiki_dir)
  
  ## Everything should have been relocated
  if( any(idx <- file.exists(file.path(vigns$dir, out))) )
    warning(paste('Files produced but not relocated:',
                  out[which(idx)]))
  
  breedR.clean_vignettes(vigns)
}

breedR.build_vignettes <- function (vigns) {
  
  
  ## Avoid duplicate names in vignettes
  dups <- duplicated(vigns$names)
  if (any(dups)) {
    names <- unique(vigns$names[dups])
    docs <- sort(basename(vigns$docs[vigns$names %in% names]))
    stop(gettextf("Detected vignette source files (%s) with shared names (%s) and therefore risking overwriting each others output files", 
                  paste(sQuote(docs), collapse = ", "), paste(sQuote(names), 
                                                              collapse = ", ")), domain = NA)
  }
  
  ## Prepare to render vignettes
  wd <- getwd()
  if (is.null(wd)) 
    stop("current working directory cannot be ascertained")
  op <- options(warn = 1)  # warnings printed as they occur
  setwd(vigns$dir)         # move to where vignettes live
  on.exit({
    setwd(wd)
    options(op)
  })  
  origfiles <- list.files(all.files = TRUE)
  file.create(".build.timestamp")
  library(rmarkdown)
  library(knitr)
  outputs <- NULL
  sourceList <- list()
  startdir <- getwd()
  for (i in seq_along(vigns$docs)) {
    file <- basename(vigns$docs[i])
    name <- vigns$names[i]
    output <- tryCatch({
      opts_chunk$set(error = FALSE)
      knit_hooks$set(purl = hook_purl)
      options(markdown.HTML.header = NULL)
      rmarkdown::render(file,
                        output_format = c('pdf_document', 'md_document'),
                        quiet = TRUE, 
                        envir = globalenv())
      
      setdiff(grep(name, list.files(), value = TRUE),
              grep(name, origfiles, value = TRUE))
    }, error = function(e) {
      stop(gettextf("processing vignette '%s' failed with diagnostics:\n%s", 
                    file, conditionMessage(e)), domain = NA, call. = FALSE)
    })
    outputs <- c(outputs, output)

  }  
  invisible(outputs)
}


breedR.move_vignettes <- function(pkg, vigns, out, wiki_dir) {
  doc_dir <- file.path(pkg$path, "inst", "doc")
  if (!file.exists(doc_dir)) {
    dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_mv <- file.path(vigns$dir,
                      grep('.*\\.(pdf|R)$', out, value = TRUE))
  out_cp <- vigns$docs
  out_wiki <- file.path(vigns$dir,
                        grep('.*(\\.md|\\.pdf|_files)$', out, value = TRUE))

  ## Updating wiki
  message("Moving ", paste(basename(out_wiki), collapse = ", "), 
          " to the Wiki")
  file.copy(out_wiki, wiki_dir, overwrite = TRUE, recursive = TRUE)
  file.copy(file.path(vigns$dir, 'img'),
            wiki_dir, overwrite = TRUE, recursive = TRUE)
  
  ## Updating packaged vignettes
  message("Moving ", paste(basename(out_mv), collapse = ", "), 
          " to inst/doc/")
  file.copy(out_mv, doc_dir, overwrite = TRUE)

  ## Updating source code in packaged vignettes
  message("Copying ", paste(basename(out_cp), collapse = ", "), 
          " to inst/doc/")
  file.copy(out_cp, doc_dir, overwrite = TRUE)
  
  ## Remove copied outcomes
  unlink(union(out_wiki, out_mv), recursive = TRUE)
  
  ## Extras
  extra_files <- devtools:::find_vignette_extras(pkg)
  if (length(extra_files) == 0) 
    return(invisible())
  message("Copying extra files ", paste(basename(extra_files), 
                                        collapse = ", "), " to inst/doc/")
  file.copy(extra_files, doc_dir, recursive = TRUE)
  invisible()
  
}

breedR.clean_vignettes <- function(vigns) {
  f <- list.files(vigns$dir, all.files = TRUE, no.. = TRUE)
  if (file.exists(".build.timestamp")) {
    newer <- file_test("-nt", f, ".build.timestamp")
    unlink(f[newer], recursive = TRUE)
    file.remove(".build.timestamp")
    return(invisible(TRUE))
  }
  invisible(FALSE)
}
