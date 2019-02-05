## Functions for the development of breedR
## Not exported, and not even packaged 
## i.e. this file is listed in .Rbuildignore

## Render vignettes
## formats: pdf for the package
##          md  for the GitHub wiki
breedR.update_guides <- function(pkg = "."){
  pkg <- as.package(pkg)
  vigns <- tools::pkgVignettes(dir = pkg$path)

  ## Filter modified vignettes
  vigns <- changed_guides(vigns)
  
  ## Check existence of any vignette
  if (length(vigns$doc) == 0) 
    return()
  
  ## Otherwise, go ahead
  message("Building ", pkg$package, " vignettes: ",
          paste(vigns$names, collapse = ", "))
  out <- breedR.build_vignettes(vigns)
  
  ## Compact vignettes
  tools::compactPDF(vigns$dir, gs_quality = "ebook")
  
  wiki_dir <- file.path(pkg$path, '..', 'breedR-wiki')
  breedR.move_vignettes(pkg, vigns, out, wiki_dir)
  
  ## Everything should have been relocated
  if( any(idx <- file.exists(file.path(vigns$dir, out))) )
    warning(paste('Files produced but not relocated:',
                  out[which(idx)]))
  
  breedR.clean_vignettes(vigns)
}


## Given a list of vignettes, filter out those who have changed
## with respect to the version in inst/doc
changed_guides <- function(vigns) {
  require(tools)
  require(utils)
  dest_dir <- file.path(vigns$pkgdir, "inst", "doc")
  
  has_changed <- function(x) {
    ## flag as changed only if x is newer than that within inst/doc
    ## and files are not identical (as per md5sum)
    dest_file <- file.path(dest_dir, basename(x))
    
    file_test("-nt", x, dest_file) &&
      !identical(unname(md5sum(x)), unname(md5sum(dest_file)))
  }
  
  subset_length <- function(x, idx) {
    if (identical(length(x), length(idx)))
      return(x[idx])
    else
      return(x)
  }
  
  change_idx <- vapply(vigns$docs, has_changed, TRUE)
  
  ans <- structure(
    lapply(vigns, subset_length, change_idx),
    class = class(vigns)
  )
  
  return(ans)
}


## This function renders md, pdf and R versions of vignettes
## Furthermore, stores images into _files directories
## It returns a vector of **new** file names (that were not there
## when the process started.) Make sure the directory is clean.
breedR.build_vignettes <- function (vigns) {
  
  
  ## Avoid duplicate names in vignettes
  dups <- duplicated(vigns$names)
  if (any(dups)) {
    names <- unique(vigns$names[dups])
    docs <- sort(basename(vigns$docs[vigns$names %in% names]))
    errmsg <- 
      sprintf("Detected vignette source files (%s) with shared names (%s)
              and therefore risking overwriting each others output files", 
              paste(sQuote(docs), collapse = ", "),
              paste(sQuote(names), collapse = ", "))
    stop(errmsg)
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
breedR.download_deps <- function(dir, type = 'all', add) {
  
  typelst <- c('source', 'win.binary', 'mac.binary', 'mac.binary.mavericks')
  type <- match.arg(type, c('all', typelst), several.ok = TRUE)
  if ('all' %in% type) type <- typelst
  
  pkgnms <- dev_package_deps(dependencies = TRUE)$package
  
  if (!missing(add)){
    pkgnms <- c(pkgnms, add)
  }
    
  for (t in type) {
    dest <- file.path(dir, t)
    dir.create(dest)
    download.packages(pkgnms, destdir = dest, type = t)
  }
  
  message('Done. ',
          'Don\'t forget to include also the ',
          'source and compiled packages for breedR.')
}


# Build graph of dependencies
# plot(breedR.deps_graph())
# fname <- '../../doc/breedR-dev/dependencies.gml'
# write_graph(breedR.deps_graph(),
#             file = fname,
#             format = 'gml')
# # in the resulting file, change 'name' by 'label'
# # in order to open it with yEd
# system(paste("sed -i 's/name/label/g'", fname))
breedR.deps_graph <- function(which = c('Depends', 'Imports')) {
  if (!require(igraph))
    stop('This requires installing igraph.')
  
  pkg <- as.package('.')
  inst <- installed.packages()
  base <- unname(inst[inst[, "Priority"] %in% c("base", "recommended"), 
                      "Package"])

  deps <- unlist(
    lapply(
      lapply(pkg[tolower(which)], parse_deps),
      `[[`, "name"),
    use.names = FALSE
  )
  
  deps <- setdiff(deps, base)
  
  # available = available.packages()
  available = inst
  
  dep_list <- list(breedR = deps)
  new_deps <- deps
  count = 0
  while (length(new_deps)>0 && count < 10) {
    all_deps <- tools::package_dependencies(new_deps, db = available)
    deps_clean <- sapply(all_deps, setdiff, base)
    deps_clean <- deps_clean[sapply(deps_clean, length) > 0]
    dep_list <- c(dep_list, deps_clean)
    new_deps <- unique(unname(unlist(deps_clean)))
    count = count + 1
  }

  dep_full <- 
    unname(
      unlist(
        sapply(names(dep_list),
               function(x) as.vector(sapply(dep_list[[x]], 
                                            function(y) c(x, y)))
        )
      )
    )
  
  return(make_graph(dep_full))
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

# Update binaries in some location (breedR web, by default)
# Only updates changed files. If it is longer than max_days
# since last download from site, then update local source.
breedR.update_bins <- function(
  source = path.expand('~/t4f/bin/PROGSF90'),
  dest   = 'web',
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

  ## Update all files from source to destinatino  
  destfn <- srcfn

  ## Complete paths
  srcfn <- file.path(source, srcfn)
  destfn <- file.path(
    switch(
      dest,
      web = path.expand('~/t4f/src/breedR-web/bin'),
      dest
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
  eqls <- sapply(tools::md5sum(destfn) == tools::md5sum(srcfn),
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

## Usage:
## repo <- file.path('file://', path.expand('~/t4f/pkgrepo'))
## available.packages(contrib.url(repo, 'source'))
## install.packages(pkgs[, 'Package'], repos = repo)
breedR.update_drat <- function(
  repodir  = '~/t4f/pkgrepo'
) {
  ## TODO: vignettes error?
  file <- devtools::build(vignettes = FALSE)
  drat::insertPackage(file, repodir)
  
  srcdir <- tempdir()
  breedR.download_deps(dir = srcdir, type = 'source')
  pkgs <- list.files(srcdir, full.names = TRUE, recursive = TRUE)
  for (fn in pkgs) {
    drat::insertPackage(fn, repodir, action = "prune")
  }
}

## Create a directory for the new release. Including a src subdir.
breedR_newrelease <- function(path = "../../breedR_releases") {
  pkg <- as.package(".")
  versiondir <- paste0("v", pkg$version)
  releasedir <- normalizePath(file.path(path, versiondir, "src"))
  if (dir.exists(releasedir)) {
    stop("Current version of breedR is already released. Bump version?")
  }
  
  dir.create(releasedir, recursive = TRUE)
  return(releasedir)
}

## Releases package versions to the web (or any other repository)
## Compares released and current versions, and releases only if necessary
## Calls build() for releasing the source and build_win() for the windows binaries
## No binary relase for mac.
breedR_release <- function(
  repodir = normalizePath('../breedR-web'),
  silent = FALSE
) {
  pkg <- as.package('.')
  pkg_basename <- paste(pkg$package, pkg$version, sep = '_')
  
  ## Check latest released versions for windows binaries and source
  ## Note: this will only look for the _installed_ R_version series
  released.win <- try(
    available.packages(contrib.url(file.path('file://', repodir), 'win.binary')),
    silent = TRUE
  )
  win.update <- 
    inherits(released.win, 'try-error') ||
    !"breedR" %in% rownames(released.win) || 
    package_version(released.win["breedR", 'Version']) < 
    package_version(pkg$version)

  released.src <- try(
    available.packages(contrib.url(file.path('file://', repodir), 'source'))
  )
  src.update <- 
    inherits(released.src, 'try-error') ||
    !"breedR" %in% rownames(released.src) || 
    package_version(released.src["breedR", 'Version']) < 
      package_version(pkg$version)
  
  ## builds and deploys packages if necessary
  if (win.update) {
    
    local_releases <- list.files(
      file.path(
        paste0(
          "../../breedR_releases/v",
          pkg$version
        )
      ),
      pattern = "zip",
      full.names = TRUE,
      recursive = TRUE
    )
    
    if (!file.exists(local_release)) {
      stop("First build a windows version with r-hub or winbuilder",
           "and put it into ../../breedR_releases")
    }
    
    for (i in seq_along(local_releases)) {
      drat::insertPackage(local_releases[[i]], repodir)
    }
    # if (!silent) message('Building windows binary pacakge ..')
    # build_win()
    # version="R-release"
    # built_path <- grep(pkg$version, list.files(contrib.url(repodir, 'source'), full.names = TRUE), value = TRUE)
    # url <- paste0("ftp://win-builder.r-project.org/", version,
    # "/", basename(built_path))
    # devtools:::upload_ftp(file = built_path, url = url)

    # if (!silent)
    #   ## Manually deploy later
    #   message('Wait for email from buildwin service and deploy with:\n',
    #           deparse(bquote(drat::insertPackage(pkg, .(repodir)))))
  } else {
    if (!silent) message('Windows binaries are up to date')
  }
  
  if (src.update) {
    if (!silent) message('Building source pacakge ..')
    src.fn <- devtools::build(args = "--compact-vignettes=\"gs+qpdf\"")
    
    src.fn <- list.files(
      file.path(
        paste0(
          "../../breedR_releases/v",
          pkg$version
        )
      ),
      pattern = "gz",
      full.names = TRUE,
      recursive = TRUE
    )
    ## deploy
    drat::insertPackage(src.fn, repodir)
  }
  
  
  if (!silent && (win.update || src.update))
    message('Finnished. Dont forget to commit breedR-web and push.')
}


## Once there are a set of compiled and source packages within
## a new version directory under reldir, deploy to repodir
## don't forget to push.
breedR_deploy <- function(
  ver = "v0.12-2",
  reldir = normalizePath("../../breedR_releases/"),
  repodir = normalizePath('../breedR-web')
) {
  dir <- file.path(reldir, ver)
  stopifnot(dir.exists(dir))
  
  subdirs <- list.dirs(dir, recursive = FALSE)
  
  for (d in subdirs) {
    fn <- grep("breedR", list.files(d), value = TRUE)
    drat::insertPackage(file.path(d, fn), repodir)
  }
}


## This function http-serves a package repository
## and displays a slide with detailed installation instructions
breedR_serve <- function(
  repo = normalizePath('~/t4f/pkgrepo'),
  pf90 = normalizePath('~/t4f/src/breedR-web/bin'),
  interface = 'wlan0'
) {
  
  get_ip <- function() {
    ipline <- system(paste('ifconfig', interface, '| grep ".*inet[^6]*:"'),
                     intern = TRUE)
    gsub("^.*?\\b(\\d{1,3}\\.\\d{1,3}\\.\\d{1,3}\\.\\d{1,3})\\b.*", "\\1", ipline)
  }
  
  ## find IP address (can I fix this?)
  http_svr <- paste0("http://", get_ip())
  pf90_svr <- file.path(http_svr, "pf90")
  
  
  ## specific instructions
  contents <- 
    c("## Install breedR and dependencies from my computer", 
      "```R",
      deparse(bquote(svr <- .(http_svr))),
      deparse(bquote(Sys.setenv(PROGSF90_URL = file.path(svr, "pf90")))),
      deparse(bquote(install.packages('breedR', repos = svr))),
      "```")
  
  fn <- tempfile('install_breedR', fileext = ".md")
  writeLines(contents, fn)
  
  ## installation instructions slide
  if (!require(rmarkdown)) stop('Install rmarkdown')
  out <- render(fn, 
                beamer_presentation(),
                output_dir = ".",
                output_file = "install_breedR.pdf")
  system(paste("evince -s", out, "&"))
  
  ## serve repo on wlan1
  # system('sudo python -m SimpleHTTPServer 80 &')
  # sin tty presente y no hay programa askpass especificado
  wd <- setwd(repo)
  on.exit(setwd(wd))
  system('python -m SimpleHTTPServer 80')
  
}

breedR_build_win <- function() {
  if (!breedR.os("windows"))
    stop("This is supposed to be run from Windows")
  if (inherits(try(breedR.bin.builtin()), 'try-error'))
    stop("Please, first run on linux:\nln -sfT binwin inst/bin")
  require(devtools)
  setwd("Z:/t4f/src/breedR")
  Sys.setenv(PROGSF90_URL = "file://z:/t4f/src/breedR-web/bin")
  build(binary=T)
}


breedR_spellcheck <- function() {
  ign_txt <- c('ai', 'AIREMLF90', 'anisotropic', 'AR', 'AUTH', 'autofill', 'backends',
               'bidimensional', 'bin', 'blupf90', 'breedr', 'BV', 'Cantet',
               'Cappa', 'Chemometrics', 'condVar', 'competition', 'cygwin',
               'dir', 'DNS', 'Dutkowski', 'edu', 'Eilers', 'em', 'Geosciences',
               'ggplots', 'Gilmour', 'Globulus', 'gpl', 'heritability', 'http',
               'Ignacy', 'indices', 'INLA', 'insim', 'IP', 'isotropic',
               'kronecker', 'Larix', 'linux', 'lme4', 'López', 'martone',
               'metagene', 'Metagene', 'microdensity', 'Misztal', 'modelframe',
               'nce', 'permanent_environmental_competition', 'nok', 'Normalise',
               'normalised', 'ord', 'org', 'os', 'param', 'passwordless',
               'pf90', 'php', 'progsf90', 'PROGSF90', 'readme', 'remlf90',
               'REMLF90', 'renderpf90', 'scp', 'sibs', 'Sinauer',
               'splineDesign', 'Trees4Future', 'uga', 'urls', 'variograms',
               'wiki', 'www')
  spell_check(pkg = '.', ignore = ign_txt)
}

# Spell checking
# (c) https://jeroenooms.github.io
spell_check <- function(pkg = ".", ignore = character()){
  pkg <- as.package(pkg)
  ignore <- c(pkg$package, hunspell::en_stats, ignore)
  
  # Check Rd manual files
  rd_files <- list.files(file.path(pkg$path, "man"), "\\.Rd$", full.names = TRUE)
  rd_lines <- lapply(sort(rd_files), spell_check_rd, ignore = ignore)
  
  # Check 'DESCRIPTION' fields
  pkg_fields <- c("title", "description")
  pkg_lines <- lapply(pkg_fields, function(x){
    spell_check_file(textConnection(pkg[[x]]), ignore = ignore)
  })
  
  # Combine
  all_sources <- c(rd_files, pkg_fields)
  all_lines <- c(rd_lines, pkg_lines)
  words_by_file <- lapply(all_lines, names)
  bad_words <- sort(unique(unlist(words_by_file)))
  
  # Find all occurences for each word
  out <- lapply(bad_words, function(word) {
    index <- which(vapply(words_by_file, `%in%`, x = word, logical(1)))
    reports <- vapply(index, function(i){
      paste0(basename(all_sources[i]), ":", all_lines[[i]][word])
    }, character(1))
  })
  structure(out, names = bad_words, class = "spellcheck")
}

print.spellcheck <- function(x, ...){
  words <- names(x)
  fmt <- paste0("%-", max(nchar(words)) + 3, "s")
  pretty_names <- sprintf(fmt, words)
  cat(sprintf(fmt, "  WORD"), "  FOUND IN\n", sep = "")
  for(i in seq_along(x)){
    cat(pretty_names[i])
    cat(paste(x[[i]], collapse = ", "))
    cat("\n")
  }
}

spell_check_text <- function(text, ignore){
  bad_words <- hunspell::hunspell(text, ignore = ignore)
  vapply(sort(unique(unlist(bad_words))), function(word) {
    line_numbers <- which(vapply(bad_words, `%in%`, x = word, logical(1)))
    paste(line_numbers, collapse = ",")
  }, character(1))
}

spell_check_file <- function(file, ignore){
  spell_check_text(readLines(file), ignore = ignore)
}

spell_check_rd <- function(rdfile, ignore){
  text <- tools::RdTextFilter(rdfile)
  spell_check_text(text, ignore = ignore)
}

download_winbuilder <- function(
  ids = c("XBOVtggM7Zax", "CL5DgvX5vRH5", "wBtX7KU9Oah3"),
  dir = normalizePath("../../breedR_releases/")
) {
  # ids <- c("XBOVtggM7Zax", "CL5DgvX5vRH5", "wBtX7KU9Oah3")
  baseurl <- "https://win-builder.r-project.org/"
  urls <- paste0(baseurl, ids, "/")

  ## pkg metadata  
  pkg <- as.package('.')
  pkg_basename <- paste(pkg$package, pkg$version, sep = '_')
  ver <- pkg$version
  
  ## create local release dir
  dir <- paste0(dir, "/v", ver)
  dir.create(dir)
  
  ## compiled package name
  fn <- paste0(pkg_basename, ".zip")
  
  for (i in seq_along(urls)) {
    
    ## get log to temporary file
    tf <- tempfile()
    download.file(file.path(urls[i], "00check.log"), tf)
    
    ## get platform
    log <- readLines(tf)
    Rver <- stringr::str_match(log[1], "R-(\\w+)")[1, 2]
    # Rver <- stringr::str_match(log[2], "R version ([\\d.]+)")[1, 2]
    winver <- stringr::str_match(log[3], "platform: ([\\w\\d_]+)")[1, 2]
    platform <- paste("windows", winver, Rver, sep = "-")

    ## create local dir    
    local_dir <- file.path(dir, platform)
    dir.create(local_dir, showWarnings = FALSE)
    
    ## get logs
    file.copy(tf, file.path(local_dir, "00check.log"))
    download.file(paste0(urls[i], "00install.out"), file.path(local_dir, "00install.out"))

    ## get built package
    remote_file <- paste0(urls[i], fn)
    download.file(remote_file, file.path(local_dir, fn))
  }
}

## Download results from r-hub
download_rhub <- function(
  x,
  dir = normalizePath("../../breedR_releases/")
) {
  x$web()   # open web logs (didn't figure out how to download them)
  baseurl <- "https://artifacts.r-hub.io/"
  ids <- x$.__enclos_env__$private$ids_
  ver <- x$.__enclos_env__$private$status_[[1]]$version
  
  urls <- paste0(baseurl, ids, "/")
  dir <- paste0(dir, "/v", ver)
  dir.create(dir)
  
  fn <- paste0("breedR_", ver, ".zip")
  
  for (i in seq_along(ids)) {
    local_dir <- file.path(dir, names(ids)[i])
    remote_file <- paste0(urls[i], fn)
    dir.create(local_dir, showWarnings = FALSE)
    try(download.file(remote_file, file.path(local_dir, fn)))
  }
}

