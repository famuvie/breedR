### Test the management of binary dependencies ###
old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

context("Binary dependencies")

# breedR.os.type() returns the current platform (broadly)
test_that(paste('breedR.os.type()', breedR.os.type(), 'identifies the platform', .Platform$OS.type), {
  check <- switch(.Platform$OS.type,
                  windows = identical(breedR.os.type(), 'windows'),
                  unix = breedR.os.type() %in% c('linux', 'mac'))
  
  expect_true(check)
})
  

# Install binaries somewhere, and check their installation
test_that('Installation of binaries and checking runs smoothly', {
  
  # If not online and test is performed in development machine
  # then set environment variable to test from local repo
  if (breedR_online() || Sys.getenv("USER") == "facu") {
    
    if (!breedR_online()){
      Sys.setenv(PROGSF90_URL = paste0("file://",
                                       normalizePath("~/t4f/src/breedR-web/bin")))
    }
    
    expect_pf90_installs <- function(os, arch) {
      tdir <- tempdir()
      path <- file.path(tdir, os)
      
      if (os == 'linux') {
        ## further specify arch for installation on linux
        path <- file.path(path, paste0(arch, 'bit'))
      }
      eval(bquote(expect_true(install_progsf90(dest = .(path),
                                               platform = .(os),
                                               arch = .(arch)))))
      
      if (os == 'windows') {
        ## further specify arch for checking on windows
        path <- file.path(path, paste0(arch, 'bit'))
      }
      eval(bquote(expect_true(check_progsf90(.(path), 
                                             platform = .(os), 
                                             quiet = TRUE))))
      
    }
    
    expect_pf90_installs('linux', '32')
    expect_pf90_installs('linux', '64')
    expect_pf90_installs('windows', '64')
    expect_pf90_installs('mac', '64')
    
  }
})

# checking somewhere else should fail
empty.dir <- tempfile()
dir.create(empty.dir)

test_that('breedR.check.bin() fails in the wrong directory', {
  expect_false(check_progsf90(empty.dir, quiet = TRUE))
})

# calling reml or aireml checks the installation
test_that('(ai)remlf90 checks the installation of binaries', {
  breedR.setOption('breedR.bin', empty.dir)
  run_dummy <- function()
    remlf90(fixed = phe_X ~ 1, data = globulus, debug = TRUE)
  
  expect_error(run_dummy(), 'Binary dependencies missing')
  
  breedR.setOption('breedR.bin', NULL)
})
