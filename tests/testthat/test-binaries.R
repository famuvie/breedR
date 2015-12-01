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
    tdir <- tempdir()
    os_arch.list <- expand.grid(os = c('windows', 'linux', 'mac'),
                                arch = c(32, 64),
                                stringsAsFactors = FALSE)
    # Windows downloads both architectures
    # Fo Mac only 64 bit available
    os_arch.list <- 
      os_arch.list[!with(os_arch.list, os != 'linux' & arch == '32'),]
    for (i in seq_len(nrow(os_arch.list))) {
      os   <- os_arch.list[i, 'os']
      arch <- os_arch.list[i, 'arch']
      
      path <- file.path(tdir, os)

      if (os == 'linux') {
        ## further specify arch for installation on linux
        path <- file.path(path, paste0(arch, 'bit'))
      }
      expect_true(install_progsf90(dest = path,
                                   platform = os,
                                   arch = arch))
      
      if (os == 'windows') {
        ## further specify arch for checking on windows
        path <- file.path(path, paste0(arch, 'bit'))
      }
      expect_true(check_progsf90(path, platform = os, quiet = TRUE))
    }
  }
  
})

# checking somewhere else should fail
empty.dir <- tempfile()
dir.create(empty.dir)
test_that('breedR.check.bin() fails in the wrong directory', {
  expect_false(check_progsf90(empty.dir, quiet = TRUE))
})


# calling reml or aireml checks the installation
breedR.setOption('breedR.bin', empty.dir)
res <- try(suppressMessages(remlf90(fixed = phe_X ~ 1,
                                    data = globulus,
                                    debug = TRUE)))
breedR.setOption('breedR.bin', NULL)

test_that('(ai)remlf90 checks the installation of binaries', {
  expect_error(res, 'Binary dependencies missing')
})
