### Test the management of binary dependencies ###
old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

context("Binary dependencies")

# breedR.os.type() returns the current platform (broadly)
check <- switch(.Platform$OS.type,
                 windows = identical(breedR.os.type(), 'windows'),
                 unix = breedR.os.type() %in% c('linux', 'mac'))

test_that(paste('breedR.os.type()', breedR.os.type(), 'identifies the platform', .Platform$OS.type), {
  expect_true(check)
})
  

# Install binaries somewhere, and check their installation
tdir <- tempdir()
os_arch.list <- expand.grid(os = c('windows', 'linux', 'mac'),
                            arch = paste0(c(32, 64), 'bit'))
os_arch.list <- 
  os_arch.list[!with(os_arch.list, os == 'mac' & arch == '32bit'),]
for(i in seq_len(nrow(os_arch.list))) {
  os   <- os_arch.list[i, 'os']
  arch <- os_arch.list[i, 'arch']
  
  path <- file.path(tdir, os, arch)

  out <- try(install_progsf90(dest = path,
                              platform = os,
                              arch = arch)
  )

  test_that(paste('Installation of binaries succeeded'), {
    expect_true(!inherits(out, 'try-error'))
  })
  
  test_that(paste('Checking installation succeeds'), {
    expect_true(check_progsf90(path, platform = os, quiet = TRUE))
  })
}

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
