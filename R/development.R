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

## Move the compiled vignettes both to breedR wiki 
## (md and pdf versions, together with associated _files)
## and inst/doc (pdf and R versions)
breedR.move_vignettes <- function(pkg, vigns, out, wiki_dir) {
  doc_dir <- file.path(pkg$path, "inst", "doc")
  if (!file.exists(doc_dir)) {
    dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  ## vignettes compiled files (to be moved)
  out_mv <- file.path(vigns$dir,
                      grep('.*\\.(pdf|R)$', out, value = TRUE))
  
  ## vignettes compiled files (to be moved to the wiki)
  out_wiki <- file.path(vigns$dir,
                        grep('.*(\\.md|\\.pdf|_files)$', out, value = TRUE))

  ## vignettes source files to be copied (not removed)
  out_cp <- vigns$docs

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
  if (file.exists(file.path(vigns$dir, ".build.timestamp"))) {
    newer <- file_test("-nt", f, ".build.timestamp")
    unlink(f[newer], recursive = TRUE)
    file.remove(file.path(vigns$dir, ".build.timestamp"))
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

# Download (recursive) dependencies
# specify where to store packages and type
breedR.download_deps <- function(dir, type) {
  
  typelst <- c('source', 'win.binary', 'mac.binary', 'mac.binary.mavericks')
  if (missing(type)) type <- getOption("pkgType")
  type <- match.arg(type, c('all', typelst), several.ok = TRUE)
  if ('all' %in% type) type <- typelst
  
  pkgnms <- dev_package_deps(dependencies = TRUE)$package
  
  for (t in type) {
    dest <- file.path(dir, t)
    dir.create(dest)
    download.packages(pkgnms, destdir = dest, type = t)
  }
}

# Downlad latest published version of binaries
# and check whether any of them have changed since last time
# Note that **it does not install** them into the package
breedR.download_bins <- function(
  baseurl  = 'http://nce.ads.uga.edu/html/projects/programs',
  basedest = '~/t4f/bin/PROGSF90',
  quiet    = FALSE
) {
  
  fnames <- unname(
    apply(
      as.matrix(
        rbind(  
          expand.grid(c('Linux'),
                      c('32bit', '64bit'),
                      c('airemlf90', 'remlf90')),
          expand.grid(c('Windows'),
                      c('32bit', '64bit'),
                      c('airemlf90.exe', 'remlf90.exe')),
          expand.grid(c('Mac_OSX'),
                      c('new'),
                      c('airemlf90', 'remlf90'))
        )
      ),
      1, function(x) do.call('file.path', as.list(x))
    )
  )

  srcf <- file.path(baseurl, fnames)
  
  destf <- file.path(
    path.expand(basedest),
    gsub('_osx/new', '', tolower(fnames)))

  # Side effect: retrieve files
  local_files <- mapply(
    retrieve_bin,
    basename(destf),
    dirname(srcf),
    dirname(destf))
  #   local_files <- path.expand(
  #     file.path(basedest,
  #               gsub('_osx/new', '', 
  #                    tolower(fnames))))
  
  ## Keep track of changes
  chglogfn <- file.path(basedest, 'changelog.rds')
  md5_new <- tools::md5sum(local_files)
  if (file.exists(chglogfn)) {
    chglog <- readRDS(chglogfn)
  } else {
    chglog <- data.frame(
      file = local_files,
      md5 = 0,
      date = 0
    )
  }
  
  ## Update hash and date for changed binaries
  if (!all(eqls <- md5_new == chglog$md5)) {
    chglog$md5 <- md5_new
    chglog$date[!eqls] <- date()
  }
  
  ## Save changelog
  ## even if nothing changed, to keep track of last check time
  saveRDS(chglog, chglogfn)

  if (!quiet) {
    if (any(!eqls))
      message(paste('Updated files:\n',
                    paste(names(!eqls), collapse = ',\n ')))
    else 
      message(paste('All binaries are up to date.'))
  }
  invisible(names(!eqls))
}

# Update binaries in some location (package bins, by default)
# Only updates changed files. If it is longer than max_days
# since last download from site, then update local source.
breedR.update_bins <- function(
  source = path.expand('~/t4f/bin/PROGSF90'),
  dest   = c('pkg', 'web'),
  max_days = 30,
  quiet = FALSE
) {
  
  dest <- match.arg(dest)
  
  ## Time since last download
  last_time <- file.mtime(file.path(source, 'changelog.rds'))
  dtime <- difftime(Sys.time(), as.Date(last_time), units = 'days')
  
  ## Download from website if necessary
  if (dtime > max_days) {
    if (!quiet) message(
      paste('It has been', floor(as.numeric(dtime, units = 'days')),
            'days since last download of binaries.\n',
            'Updating local versions...')
    )
    breedR.download_bins(basedest = source)
  }
  
  ## Last local version of binaries
  srcfn <- list.files(source, recursive = TRUE)
  srcfn <- srcfn[-grep('changelog\\.rds', srcfn)]
  
  ## The package ships only 32bit versions
  if (dest == 'pkg') {
    idx <- !grepl('64bit', srcfn)
    srcfn <- srcfn[idx]
    destfn <- gsub('/32bit', '', srcfn)
  } else {
    destfn <- srcfn
  }

  ## Complete paths
  srcfn <- file.path(source, srcfn)
  destfn <- file.path(
    switch(
      dest,
      pkg = system.file('bin', package = 'breedR'),
      web = path.expand('~/t4f/src/breedR-web/bin')
    ),
    destfn
  )
  
  if( any(dir <- !dir.exists(destd <- unique(dirname(destfn)))) ) {
    res <- sapply(destd[dir], dir.create, recursive = TRUE)
    if (!all(res))
      stop(paste('Had a problem while creating\n',
                 paste(destd[dir], collapse = ',\n ')))
    
  }
  
  ## New versions of binaries
  ## isTRUE gives FALSE for NAs
  eqls <- sapply(tools::md5sum(srcfn) == tools::md5sum(destfn),
                 isTRUE)

  if (!quiet) {
    if (any(!eqls)) {
      message(paste('Updated files:\n',
                    paste(names(!eqls), collapse = ',\n ')))
      res <- file.copy(srcfn[!eqls], destfn[!eqls], overwrite = TRUE)
      if (!isTRUE(all(res))) {
        stop(paste('Error copying files',
                   paste(names(!eqls)[is.na(res)], collapse = ',\n ')))
      }
    } else 
      message(paste('All binaries are up to date.'))
  }
  
  invisible(names(!eqls))
}