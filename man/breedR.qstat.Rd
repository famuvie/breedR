\name{breedR.qstat}
\alias{breedR.qdel}
\alias{breedR.qget}
\alias{breedR.qnuke}
\alias{breedR.qstat}
\alias{print.breedR.q}
\alias{summary.breedR.q}
\title{Control and view a remote breedR-queue of submitted jobs}
\usage{
summary.breedR.q(object, ...)

breedR.qget(id, remove = TRUE)

breedR.qdel(id, statlst)

breedR.qstat(id)

breedR.qnuke()
}
\arguments{
  \item{id}{The job-id which is the output from
  \code{breedR} when the job is submitted, the job-number
  or job-name. For \code{breedR.qstat}, \code{id} is
  optional and if omitted all the jobs will be listed.}

  \item{remove}{Logical. If FALSE, leave the job on the
  server after retrival, otherwise remove it (default).}

  \item{x,}{object An \code{breedR.q}-object which is the
  output from \code{breedR.qstat}}

  \item{...}{other arguments

  \code{breedR.qstat} show job(s) on the server,
  \code{breedR.qget} fetch the results (and by default
  remove the files on the server), \code{breedR.qdel}
  removes a job on the server and \code{breedR.qnuke}
  remove all jobs on the server.

  The recommended procedure is to use \code{r=breedR(...,
  breedR.call="submit")} and then do
  \code{r=breedR.qget(r)} at a later stage. If the job is
  not finished, then \code{r} will not be overwritten and
  this step can be repeated.  The reason for this
  procedure, is that some information usually stored in the
  result object does not go through the remote server,
  hence have to be appended to the results that are
  retrieved from the server. Hence doing
  \code{r=breedR(..., breedR.call="submit")} and then later
  retrive it using \code{r=breedR.qget(1)}, say, then
  \code{r} does not contain all the usual information.  All
  the main results are there, but administrative
  information which is required to call
  \code{breedR.hyperpar} or \code{breedR.rerun} are not
  there.}
}
\value{
\code{breedR.qstat} returns an \code{breedR.q}-object with
information about current jobs.
}
\description{
Control and view a remote breedR-queue of submitted jobs
}
\examples{
\dontrun{
r = remlf90(y~1, data = data.frame(y=rnorm(10)), breedR.call = "submit")
summary(r)   # same as breedR.qstat(r)
breedR.qstat()
r = breedR.qget(r, remove=FALSE)
breedR.qdel(1)
breedR.qnuke()
summary(r)   # results of the analysis
}
}
\seealso{
\code{\link{remlf90}}
}

