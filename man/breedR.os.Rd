\name{breedR.os}
\alias{breedR.os}
\title{Check host OS}
\usage{
breedR.os(type = c("linux", "mac", "windows", "else"))
}
\description{
Identifies host operating system.
}
\details{
Relies on \code{.Platform$OS.type, but distinguishes
between linux or mac.
}

