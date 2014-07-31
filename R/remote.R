# Functions adapted from the R-INLA project
# www.r-inla.org

#' Control and view a remote breedR-queue of submitted jobs
#' 
#' @name breedR.qstat
#' @aliases breedR.qstat breedR.qget breedR.qdel breedR.qnuke summary.breedR.q 
#'   print.breedR.q
#' @export breedR.qstat breedR.qget breedR.qdel breedR.qnuke summary.breedR.q 
#'   print.breedR.q
#'   
#'   
#' @param id The job-id which is the output from \code{breedR} when the job is 
#'   submitted,  the job-number or job-name. For \code{breedR.qstat}, \code{id} 
#'   is optional and if omitted all the jobs will be listed.
#' @param remove Logical. If FALSE, leave the job on the server after retrival, 
#'   otherwise remove it (default).
#' @param x, object An \code{breedR.q}-object which is the output from 
#'   \code{breedR.qstat}
#' @param ... other arguments
#'   
#'   \code{breedR.qstat} show job(s) on the server, \code{breedR.qget} fetch the
#'   results (and by default remove the files on the server), 
#'   \code{breedR.qdel} removes a job on the server and \code{breedR.qnuke}
#'   remove all jobs on the server.
#'   
#'   The recommended procedure is to use \code{r=breedR(..., 
#'   breedR.call="submit")} and then do \code{r=breedR.qget(r)} at a later
#'   stage. If the job is not finished, then \code{r} will not be overwritten
#'   and this step can be repeated.  The reason for this procedure, is that some
#'   information usually stored in the result object does not go through the 
#'   remote server, hence have to be appended to the results that are retrieved 
#'   from the server. Hence doing \code{r=breedR(..., breedR.call="submit")} and
#'   then later retrive it using \code{r=breedR.qget(1)}, say, then \code{r}
#'   does not contain all the usual information.  All the main results are
#'   there, but administrative information which is required to call
#'   \code{breedR.hyperpar} or \code{breedR.rerun} are not there.
#'   
#' @return \code{breedR.qstat} returns an \code{breedR.q}-object with
#'   information about current jobs.
#'   
#' @seealso \code{\link{remlf90}}
#' @examples
#' \dontrun{
#' r = remlf90(y~1, data = data.frame(y=rnorm(10)), breedR.call = "submit")
#' summary(r)   # same as breedR.qstat(r)
#' breedR.qstat()
#' r = breedR.qget(r, remove=FALSE)
#' breedR.qdel(1)
#' breedR.qnuke()
#' summary(r)   # results of the analysis
#' }

#' @rdname breedR.qstat
`summary.breedR.q` = function(object, ...)
{
  print(object, ...)
}

`print.breedR.q` = function(x, ...)
{
  if (length(x) == 0) {
    ##cat("No jobs available\n")
  } else {
    for(k in seq_along(x)) {
      cat("\t Job:", x[[k]]$no, "\tId:", x[[k]]$id, "\tStatus:", x[[k]]$status, "\n")
    }
  }
  return (invisible(x))
}

#' @rdname breedR.qstat
`breedR.qget` = function(id, remove = TRUE)
{

  stopifnot( !missing(id) )
  stopifnot( inherits(id, 'remlf90') & is.list(id) & exists('id', id) )
  
  # check that the job is correctly finished and uniquely determined
  statlst <- breedR.qstat(id)
  if( length(statlst) == 0 ) stop('Job not found')
  if( length(statlst) != 1 ) stop('This should not happen')
  status <- statlst[[1]]
  if( id$id != status$id ) stop('This should not happen')
  if( status$status != "Finished" ) {
    print(statlst)
    stop('Job not finished')
  }
  
  # Remote target dir
  rdir = file.path('tmp', '.breedR.remote',
                   paste('breedR-job-', status$id, sep = ''))
  
  ldir <- retrieve_remote(rdir)

  # Integrate the model structure with the results
  ans <- parse_results(file.path(ldir, 'solutions'),
                       id$effects,
                       id$mf,
                       readLines(file.path(ldir, 'LOG')),
                       id$method,
                       id$mcout)
  class(ans) <- c('breedR', 'remlf90')  
  
  if( remove ) suppressMessages(breedR.qdel(id, statlst))
  
  message('Job retrieved')
  return (ans)
}

#' @rdname breedR.qstat
`breedR.qdel` = function(id, statlst)
{
  
  if( missing(id) ) stop('No job specified. To delete all jobs use breedR.qnuke()')

  if( missing(statlst) ) statlst <- breedR.qstat(id)
  if( length(statlst) == 0 ) stop('Job not found')
  if( length(statlst) != 1 ) stop('This should not happen')
  status <- statlst[[1]]
  
  
  # Remote target dir
  rdir = file.path('tmp', '.breedR.remote',
                   paste('breedR-job-', status$id, sep = ''))
  
  ssh_commands <- paste('rm -rf', rdir)
  
  # If process is running, then kill him
  if( grepl("Running", status$status) & status$pid != 0) {
    ssh_commands <- c(ssh_commands,
                      paste('kill', status$pid))
  } else {
    if( status$pid != 0 ) stop('This should not happen')
  }
  
  # Execute
  res <- breedR.ssh(ssh_commands, intern = TRUE)
  if( !is.character(res) & length(res) != 0) stop('This should not happen')
  
  message(paste('Deleted job:', status$no, 'Id:', status$id))
}

#' @rdname breedR.qstat
`breedR.qstat` = function(id)
{
  # id can be a list of class remlf90 widh an item $id with a job id
  # or a job number (integer > 0)
  # or either NULL, missing or negative, in which case everything will be listed
  
  rdir='tmp/.breedR.remote'
  
  if( !missing(id) ) {
    if( is.list(id) ){
      if( exists('id', id) ) {
        id <- id$id                # Keep only the id
      } else stop('No id found')
    }
    if( is.numeric(id) ) idx = as.integer(id)  # List only job number idx
    else idx = -1  # Never match job number. Try with id.
  } else {
    idx = 0 # List all jobs
    id  = "NULL"
  }
  
  ssh_commands <- 
    c(paste('mkdir -p', rdir),      # Make sure the temp dir exists
      paste('cd', rdir),            # Move in
      'nno=0',                      # Reset job counter
      'for d in \\$(ls -1 .)',      # For each directory in the temporary dir
      'do if [ -d \\$d -a \\! -f \\$d/working ]', # Sanity check
      'then nno=\\$[ \\$nno + 1 ]', # Increase job counter
      paste('if [', idx, '-eq 0 -o',# Check whether the current job
            idx, '-eq \\$nno -o',   #    should be listed
            '\\$d =', paste('breedR-job-', id, sep =''), ']'),
      'then myid=\\$(echo \\$d | sed \'s/breedR-job-//\')', # get current job id,
      'mypid=0; myppid=0',                    # reset PID and PPID numbers
      'if [ -f \\$d/done ]',   # If there exists a file named done
      'then status="Finished"',     #   the job is finished
      'else myppid=\\$(ps -o pid,command -C bash | grep -s \\$myid | grep -v grep | awk \'{print int(\\$1)}\')',  # find PPID of script, excluding the match of this grep process itself
      'if [ \\$myppid ]',            # If there is such PID
      'then if [ \\$(ps --ppid \\$myppid -o comm= | grep -s remlf90) ]', # If reml is running under that PID
      'then mypid=\\$(ps --ppid \\$myppid -o pid=)',
      'runtime=\\$(ps --ppid \\$myppid -o time=)',
      'status="\\"Running\\(\\$runtime\\)\\""',   # Report running time
      'else status="Aborted"',      # Otherwise, it must have failed
      'fi',                         # End reml running
      'else myppid=0',
      'status="Aborted"',      # No active PID but job not done !!
      'fi',                         # End there is such pid
      'fi',                         # End job is done, running or aborted?
      'echo "\\$myid \\$nno \\$mypid \\$myppid \\$status"', # Result
      'fi',                         # End case the current job should be listed
      'fi',                         # End sanitized case
      'done')                       # End For

  # execute ssh script
  out <- breedR.ssh(ssh_commands, intern = TRUE)

  # Parse results
  if (length(out) >= 1 && nchar(out[1]) > 0) {
    output = lapply(strsplit(out, " +"),
                    function(a) {
                      names(a) = c("id", "no", "pid", "ppid", "status")
                      return(as.list(a))
                      })
    # Diagnostic checks
    runners.idx <- which(sapply(output, function(x) grepl("Running", x$status)))
    actives.idx <- which(sapply(output, function(x) x$pid != "0"))
    
    # There should not be any non-running process with PID>0
    # nor any running process with PID=0
    if( !identical(runners.idx, actives.idx) ) {
      review.idx <- sort(union(setdiff(runners.idx, actives.idx),
                               setdiff(actives.idx, runners.idx)))
      warning('Either some running processes have PID = 0 or some finish/aborted process still active. Please check.')
    }
    
  } else {
    output = list()
  }
  class(output) = "breedR.q"
  
  
  return (output)
}

#' @rdname breedR.qstat
`breedR.qnuke` = function()
{
  # Remote dir
  rdir = file.path('tmp', '.breedR.remote')
  
  ssh_commands <- c(paste('rm -rf', rdir))
  
  # Processes to be killed
  statlst <- breedR.qstat()
  runners.idx <- which(sapply(statlst, function(x) grepl("Running", x$status)))
  
  # We kill all the runners with pid > 0 (which should be all)
  killpids <- sapply(statlst[runners.idx], function(x) as.integer(x$pid))
  killpids <- killpids[killpids > 0L]   # Just in case
  killpids <- paste(killpids, collapse = ' ')
    
  if( length(killpids) > 0 ) {
    ssh_commands <- c(ssh_commands,
                      paste('kill', killpids))
  } 
  
  # Execute
  res <- breedR.ssh(ssh_commands, intern = TRUE)
  if( !is.character(res) & length(res) != 0) stop('This should not happen')
  
  message('NUKE')
}



### Ancillary functions ###
### non-exported        ###


#' Retrieve ssh configuration parameters
breedR.ssh_params <- function(format = c('string', 'list')) {
  
  format <- match.arg(format)
  
  ssh <- breedR.getOption(c('remote.host',
                            'remote.user',
                            'remote.port',
                            'ssh.options'))
  
  # Options for ssh
  ssh_params <- paste(paste('-p', ssh$remote.port, sep = ''),
                      ssh$ssh.options,
                      paste(ssh$remote.user, '@', ssh$remote.host, sep = ''))
  
  if( format == 'string' )
    return(ssh_params)
  else
    return(ssh)
}

#' Perform an SSH system call
#' 
#' Use the given connection parameters for and run the given commands remotely. 
#' It also admits some pre or post strings (pipelines, or other modifications)
breedR.ssh <- function(commands,
                       params = breedR.ssh_params(),
                       pre,
                       post,
                       ...) {
  cmd_str <- paste(commands, collapse = '; ')
  call_str <- paste('ssh', params, '"', cmd_str, '"')
  if( !missing(pre) ) {
    call_str <- paste(pre, call_str)
  }
  if( !missing(post) ) {
    call_str <- paste(call_str, post)
  }
  
  system(call_str, ...)
}


#' Perform a job remotely
breedR.remote = function(jobid, breedR.call, verbose = TRUE)
{
  if( verbose ) {
    message(paste('Run', breedR.call, 'at host',
                  breedR.ssh_params('list')$remote.host))
  }
  
  # Remote directory
  rdir=file.path('tmp',
                 '.breedR.remote',
                 paste('breedR-remote', jobid, sep = '-'))
  
  
  # To be executed on the server
  ssh_commands <- c(paste('mkdir -p', rdir),        # make temp dir for job
                    paste('cd', rdir),              # switch to job dir
                    'tar xfmz -',                   # uncompress stuff
                    'echo parameters > interface',  # interface arguments
                    paste(breedR.call,
                          '< interface',
                          '> LOG 2>&1') # run breedR
  )
  
  # Compress stuff and execute ssh commands
  res <- breedR.ssh(ssh_commands,
                    pre = 'tar cfmz - . |')
  
  if( verbose ) {
    message(paste(' *** Computations finished at', date(),
                  '\n *** Transfer the results...'))
  }
  
  Sys.sleep(1)  # Not too fast...
  # Retrieve results to local dir
  ldir <- retrieve_remote(rdir)
  
  return(ldir)
}


#' Submit a job with a given id and remote program call
breedR.submit <- function(jobid, breedR.call) {
  
  # Remote directory
  rdir=file.path('tmp',
                 '.breedR.remote',
                 paste('breedR-job', jobid, sep = '-'))
  
  
  # To be executed on the server
  ssh_commands <- c(paste('mkdir -p', rdir),        # make temp dir for job
                    paste('cd', rdir),              # switch to job dir
                    'touch working',                # flag file
                    'tar xfmz -',                   # uncompress stuff
                    'echo parameters > interface',  # interface arguments
                    paste('{ rm -f working; ',      # run reml
                          breedR.call,
                          '< interface && touch done; }',
                          '</dev/null > LOG 2>&1 &'))
  
  # Compress local stuff and execute ssh commands *in background*
  breedR.ssh(ssh_commands,
             pre = 'tar cfmz - . |',
             post = '&')
}



#' Retrieve results stored in some remote directory
#' 
#' Use scp to transfer compressed files. Clean up afterwards.
#' @return dir name where the results are retrieved
retrieve_remote <- function (rdir) {
  # Compressed filename for storing results remotely
  tarfile = tempfile(pattern = 'results',
                     tmpdir = '..',
                     fileext = '.tar')
  
  # Save results into the compressed file
  ssh_commands <- 
    c(paste('cd', rdir),           # Move in
      paste('tar cf',              # Compress results into a tar file
            tarfile,
            'LOG solutions'))
  res <- breedR.ssh(ssh_commands, intern = TRUE)
  if( !is.character(res) & length(res) != 0) stop('This should not happen')
  
  # Copy the compressed file to local
  tf <- tempfile(pattern = 'breedR.result_', fileext = '.tar')
  ssh_par <- breedR.ssh_params('list')
  scp_args <- paste('-P', ssh_par$remote.port, ' -B -C -p -q', sep ='')
  scp_file <- paste(ssh_par$remote.user, '@', ssh_par$remote.host, ':',
                    file.path(dirname(rdir), basename(tarfile)), sep = '')
  res <- system(paste('scp', scp_args, scp_file, tf))
  stopifnot( identical(res, 0L) )
  
  system(paste('tar xfm', tf, '-C', dirname(tf))) # Uncompress
  unlink(tf)                                      # Remove tar
  
  # Cleanup remote temporary tar
  ssh_commands <- paste('rm', file.path(dirname(rdir), basename(tarfile)))
  res <- breedR.ssh(ssh_commands, intern = TRUE)
  if( !is.character(res) & length(res) != 0) stop('This should not happen')
  return(dirname(tf))
}
