### Test the management of binary dependencies ###

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
os.list <- c('windows', 'linux', 'mac')
for(os in os.list) {
  path <- file.path(tdir, os)

  out <- breedR.install.bin(path = path,
                            os.type = os)

  test_that(paste('Installation of binaries succeeded'), {
    expect_true(all(out))
  })
  
  test_that(paste('Checking installation succeeds'), {
    expect_true(breedR.check.bin(path))
  })
}

# checking somewhere else should fail
empty.dir <- tempfile()
dir.create(empty.dir)
test_that('breedR.check.bin() fails in the wrong directory', {
  expect_false(breedR.check.bin(empty.dir, silent = TRUE))
})


# calling reml or aireml checks the installation

breedR.setOption('breedR.bin', empty.dir)
res <- try(suppressWarnings(remlf90(fixed = phe_X ~ 1,
                                    data = globulus,
                                    debug = TRUE)),
           silent = TRUE)
breedR.setOption('breedR.bin', NULL)

test_that('(ai)remlf90 checks the installation of binaries', {
  expect_true(inherits(res, 'try-error'))
})
