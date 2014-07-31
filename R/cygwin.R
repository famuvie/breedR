## Nothing to export

`breedR.cygwin.check.path` = function(path = breedR.getOption("cygwin"))
{
    return (file.exists(path) && file.info(path)$isdir)
}

`breedR.cygwin.run.command` = function(command, file.log = NULL, ...)
{
    if (breedR.cygwin.check.path()) {
        cmd = paste(breedR.getOption("cygwin"),
                     "/bin/bash.exe -c ",
                     shQuote(paste("export PATH=/bin:/usr/bin:$PATH;",
                                   command,
                                   breedR.ifelse(is.null(file.log), "", paste(" > ", file.log)), 
                                   sep=""), type="cmd"), sep="")
        system(cmd, ...)
    } else {
        stop(paste("Cannot find the CYGWIN installation:", breedR.getOption("cygwin")))
    }
}

`breedR.cygwin.map.filename` = function(filename, windows2cygwin = TRUE)
{
    if (windows2cygwin)
        return (paste(breedR.cygwin.run.command(paste("cygpath -u ", filename), intern = TRUE, ignore.stderr = TRUE),
                      collapse = ' '))
    else
        return (paste(breedR.cygwin.run.command(paste("cygpath -w -m -s ", filename), intern = TRUE, ignore.stderr = TRUE),
                      collapse= ' '))
}


        
