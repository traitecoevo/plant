##' Load reference output from the previous version of the model.
##'
##' This loads the files created by the Falster version of the EBT
##' model into R.  There are a lot of different files, one per
##' variable.  Eventually, we will have functions for working with the
##' output and getting the things ready for comparison.
##'
##' @title Load Reference Output
##' @param path Path (relative or absolute) to directory containing output.
##' @return A list with elements (undocumented at present).
##' @author Rich FitzJohn
##' @keywords internal
load.reference.output <- function(path) {
  ## First, determine which species/strategy indicies are present.
  ## These will be in files named "<idx>_age.txt"
  re <- "^([0-9]+)_age.txt$"
  strategy.idx <- as.integer(sub(re, "\\1", dir(path, pattern=re)))

  ## Files for which we will build a matrix.  The actual filename is
  ## <idx>_<type>.txt, with type from this vector:
  files.mat.spp <- c("age", "coh_m", "coh_n", "coh_r",
                     "bound_m", "bound_s", "bound_n", "bound_r",
                     "d_bound_m", "d_bound_s", "d_bound_r")
  ## And files for which we will build a data.frame
  files.df.spp <- "popn"
  files.spp <- c(files.mat.spp, files.df.spp)

  ## Load a single strategy object, with index 'idx':
  load.strategy <- function(idx) {
    fmt <- sprintf("%d_%%s.txt", idx)
    ret <- lapply(files.spp, function(x)
                  read.tsv.matrix(file.path(path, sprintf(fmt, x)),
                                  x %in% files.df.spp))
    names(ret) <- files.spp
    ret
  }

  ## Load the simulation parameters (except for the four key traits).
  params <- sub(";$", "", readLines(file.path(path, "params.m"))[-1])
  params <- do.call(rbind, strsplit(params, "="))
  params <- structure(as.list(params[,2]),
                      names=gsub("^p\\.", "", params[,1]))

  ## Load global traits:
  stand <- read.tsv.matrix(file.path(path, "stand.txt"), TRUE)
  fitness <- c(read.tsv.matrix(file.path(path, "fitness.txt"), FALSE))
  patch_age <- read.tsv.matrix(file.path(path, "patch_age.txt"), TRUE)

  list(stand=stand,
       fitness=fitness,
       patch_age=patch_age,
       params=params,
       popn=lapply(strategy.idx, load.strategy))
}

## Read a file containing tab-separated values as a matrix.
##
## This is substantially faster than \code{read.table} because no type
## checking is done.  However, it also includes far fewer tests (e.g.,
## it will fail badly on 1 line matrices).
read.tsv.matrix <- function(filename, with.header=FALSE) {
  con <- file(filename, "rt")
  on.exit(close(con))
  if (with.header)
    h <- scan(con, nlines=1, quiet=TRUE, what=character())

  d1 <- scan(con, nlines=1, quiet=TRUE, sep="")
  d2 <- scan(con, quiet=TRUE, sep="")
  m <- matrix(c(d1, d2), ncol=length(d1), byrow=TRUE)

  if (with.header)
    colnames(m) <- h

  m[is.nan(m)] <- NA
  m
}
