## ----setup, echo = FALSE, include=FALSE----------------------------------
library(breedR)

## ----check-pedigree------------------------------------------------------
set.seed(123); n.ped <- 5
ped.nightmare <- matrix(sample(30, n.ped*3), n.ped, 3,
                        dimnames = list(NULL, c('self', 'sire', 'dam')))
check_pedigree(ped.nightmare)

## ----build-pedigree------------------------------------------------------
ped.fix <- build_pedigree(1:3, data = ped.nightmare)
check_pedigree(ped.fix)
attr(ped.fix, 'map')  # map from old to new codes

## ----compare-pedigrees, results = 'asis', echo=FALSE---------------------
knitr::kable(ped.nightmare)
knitr::kable(as.data.frame(ped.fix))

## ----exercise-pedigree-dat-----------------------------------------------
test.dat <- data.frame(ped.nightmare, y = rnorm(n.ped))
res.raw <- remlf90(fixed   = y ~ 1,
                   genetic = list(model = 'add_animal',
                                  pedigree = ped.nightmare,
                                # pedigree = test.dat[, 1:3],  # same thing
                                  var.ini = 1,
                                  id = 'self'),
                   var.ini = list(resid = 1),
                   data    = test.dat)

## pedigree has been recoded!
length(ranef(res.raw)$genetic)

## The pedigree used in the model matches the one manually built
identical(ped.fix, get_pedigree(res.raw))

## ----PBV-----------------------------------------------------------------
## Predicted Breeding Valuess of the observed individuals
## Left-multiplying the vector of BLUP by the incidence matrix
## gives the BLUP of the observations in the right order.
Za <- model.matrix(res.raw)$genetic  # incidence matrix
gen.blup <- with(ranef(res.raw),
                 cbind(value=genetic,
                       's.e.'=attr(genetic, 'se')))
PBVs <- Za %*% gen.blup
rownames(PBVs) <- test.dat$self

## ----PBV-table, results = 'asis', echo=FALSE-----------------------------
knitr::kable(PBVs, digits = 2)

## ----PBV-founders--------------------------------------------------------
## original codes of non-observed parents
(founders.orig <- setdiff(
  sort(unique(as.vector(ped.nightmare[, c("sire", "dam")]))),
  ped.nightmare[, "self"]
))

## map from original to internal codes
map.codes <- attr(get_pedigree(res.raw), "map")

## internal codes of non-observed parents
founders.int <- map.codes[founders.orig]

## Breeding Values of non-observed parents
founders.PBVs <- gen.blup[founders.int, ]
rownames(founders.PBVs) <- founders.orig

## ----PBV-founders-table, results = 'asis', echo=FALSE--------------------
knitr::kable(founders.PBVs, digits = 2)

## ----reverse-lookup------------------------------------------------------
## individuals of interest in internal codification
idx <- c(3, 5, 9)

## original codes
(match(idx, map.codes))

