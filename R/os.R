# Functions borrowed from the R-INLA project
# www.r-inla.org

## nothing to export

#' Check host OS
#' 
#' Identifies host operating system.
#' 
#' Relies on \code{.Platform$OS.type, but distinguishes between linux or mac.
`breedR.os` = function(type = c("linux", "mac", "windows", "else"))
{
    if (missing(type)) {
        stop("Type of OS is required.")
    }
    type = match.arg(type)
    
    if (type == "windows") {
        return (.Platform$OS.type == "windows")
    } else if (type == "mac") {
        result = (file.info("/Library")$isdir && file.info("/Applications")$isdir)
        if (is.na(result)) {
            result = FALSE
        }
        return (result)
    } else if (type == "linux") {
        return ((.Platform$OS.type == "unix") && !breedR.os("mac"))
    } else if (type == "else") {
        return (TRUE)
    } else {
        stop("This shouldn't happen.")
    }
}
`breedR.os.type` = function()
{
    for (os in c("windows", "mac", "linux", "else")) {
        if (breedR.os(os)) {
            return (os)
        }
    }
    stop("This shouldn't happen.")
}

## this seems be to a nice way to test 32/64 bits architecture.
`breedR.os.32or64bit` = function()
{
    return (ifelse(.Machine$sizeof.pointer == 4, "32", "64"))
}
`breedR.os.is.32bit` = function()
{
    return (breedR.os.32or64bit() == "32")
}
`breedR.os.is.64bit` = function()
{
    return (breedR.os.32or64bit() == "64")
}
