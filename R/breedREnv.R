## Nothing to export

## Define the environment for breedR used to store options and the
## model-list. Reuse the environment if it is there already.

## Thanks to the INLA team, from where I took the whole option management system.

if (exists(".breedREnv") && is.environment(.breedREnv)) {
    ## then reuse it
} else {
    .breedREnv = new.env()
}

`breedR.get.breedREnv` = function(...)
{
    if (exists(".breedREnv") && is.environment(.breedREnv))
        return (.breedREnv)
    stop("Environment '.breedREnv' does not exists and is required for breedR to work. Restart 'R'.")
}

