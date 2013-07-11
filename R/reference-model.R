##' Load and save parameters from the reference model
##' (falster-traitdiversity).
##'
##' These functions allow loading of the reference data set.  They are
##' known to work with version d649534 of falster-traitdiversity (and
##' will not work with any earlier version).
##'
##' The main constraint that must be satisfied on the tree->reference
##' implementation is that \emph{only} lma and hmat may vary.  All
##' other strategy will be checked to ensure that they are identical
##' across strategies.
##'
##' The reference implementation expects the core traits on a log
##' scale, and also `mean_disturbance_interval`.  We do this
##' translation immediately on reading and writing the objects, so
##' that all traits are on a natural scale within the R code.
##'
##' These functions are inverses of each other:
##' \code{save.reference.parameters(load.reference.parameters(path),
##' path)} should leave the directory \code{path} unchanged, aside
##' from possible changes in precision due to read/writing numbers.
##'
##' @title Load and Save Reference Model Parameters
##' @param path The path (absolute, or relative to the working
##' directory) to the directory containing files "Stoch.txt" and
##' "params.m" (or, in the case of \code{save.reference parameters},
##' the directory that \emph{will} contain those files).
##' @return \code{load.reference.parameters} returns
##' \code{\link{Parameters}} object.
##' @rdname load.reference.parameters
##' @seealso \code{\link{load.reference.output}}, which loads the
##' result of running falster-traitdiversity on a set of parameters.
##' @export
##' @author Rich FitzJohn
load.reference.parameters <- function(path)
  parameters.from.reference(read.reference.parameters(path))
##' @export
##' @rdname load.reference.parameters
##' @param p A \code{\link{Parameters}} object.
##' @return \code{save.reference.parameters} returns nothing - called
##' for it's side effect.
save.reference.parameters <- function(p, path)
  write.reference.parameters(reference.from.parameters(p), path)

##' Load the output of running the reference model
##' (falster-traitdiversity).
##'
##' There is still work to do in sanitising the output ready for
##' comparison, so the return value will change.
##'
##' This function is known to work with version d649534 of
##' falster-traitdiversity (and will not work with any earlier
##' version).
##'
##' @title Load Reference Model Output
##' @param path The path containing the output; this is the directory
##' that contains the files "Stoch.txt", "params.m", and the directory
##' "res" that has the actual output files.
##' @return A list with components:
##' \item{stand}{Not sure!}
##' \item{patch_age}{Not sure!}
##' \item{seed_rain_out}{Vector of output seed rain, one entry per
##' strategy/species}
##' \item{popn}{A list of lists, one per strategy.  Each strategy's
##' list has components:
##' \describe{
##'   \item{age}{Not sure!}
##'   \item{coh_m}{Not sure!}
##'   \item{coh_n}{Not sure!}
##'   \item{coh_r}{Not sure!}
##'   \item{bound_m}{Not sure!}
##'   \item{bound_s}{Not sure!}
##'   \item{bound_n}{Not sure!}
##'   \item{bound_r}{Not sure!}
##'   \item{d_bound_m}{Not sure!}
##'   \item{d_bound_s}{Not sure!}
##'   \item{d_bound_r}{Not sure!}
##'   \item{popn}{Not sure!}
##' }}
##' \item{parameters}{A \code{\link{Parameters}} object}
##' @author Rich FitzJohn
##' @export
load.reference.output <- function(path) {
  p <- load.reference.parameters(path)
  check.reference.output(p$size, path)

  path.output <- file.path(path, "res")

  ## Load population traits:
  popn.files <- c("stand.txt", "patch_age.txt", "Strategy.txt")
  popn <- lapply(file.path(path.output, popn.files), read.tsv.matrix,
                   TRUE)
  names(popn) <- tolower(sub("\\.txt", "", popn.files))

  ## Check that the output is plausible given input parameters:
  traits.input <- t(sapply(p$strategies, function(x)
                           unlist(x$parameters[reference.traits.names()])))
  traits.output <- popn$strategy[, sprintf("Tr%d", 0:3)]
  ok <- isTRUE(all.equal(unname(traits.input), unname(traits.output),
                         tolerance=1e-5))
  if (!ok)
    stop("Input and output traits appear to vary")
  if (!isTRUE(all.equal(popn$strategy[, "X_start"], p$seed_rain)))
    stop("Input and output seed_rain appear to vary")

  popn$seed_rain_out <- popn$strategy[, "X_end"]
  popn <- popn[names(popn) != "strategy"]

  ## Load a single strategy object, with index 'idx':
  load.strategy <- function(idx) {
    fmt <- sprintf("%d_%%s.txt", idx)
    ret <- lapply(reference.output.files(), function(x)
                  read.tsv.matrix(file.path(path.output, sprintf(fmt, x)),
                                  x == "popn"))
    names(ret) <- reference.output.files()
    ret
  }
  popn$strategies <- lapply(seq(0, length.out=p$size), load.strategy)
  popn$parameters <- p

  popn
}

##' Run the reference model, falster-traitdiversity
##'
##' This function is known to work with version d649534 of
##' falster-traitdiversity (and will not work with any earlier
##' version).
##'
##' @title Run Reference Model
##' @param path Path to directory containing parameter files (see
##' \code{\link{load.reference.output}} and
##' \code{\link{load.reference.parameters}}
##' @param p Optional \code{\link{Parameters}} object that will be
##' written to \code{path} before running model (equivalent to "run
##' model on these parameters in this directory")
##' @param path.evolve Path to the "evolve" program.  If not
##' specified, "evolve" is assumed to the in the \code{PATH}.
##' @param verbose Logical indicating if the model output should be
##' printed (TRUE by default).
##' @author Rich FitzJohn
##' @export
run.reference <- function(path, p=NULL, path.evolve=NULL, verbose=TRUE) {
  if (is.null(path.evolve)) {
    program <- "evolve"
  } else {
    if (length(path.evolve) != 1)
      stop("Expected scalar path to the evolve program")
    program <- file.path(normalizePath(path.evolve), "evolve")
  }

  if (is.null(p)) {
    files <- c("params.m", "Stoch.txt")
    if (!all(file.exists(file.path(path, files))))
      stop(sprintf("Could not find file(s) in %s: %s", path,
                   paste(files, collapse=", ")))
  } else {
    save.reference.parameters(p, path)
  }

  ## We need to change to the directory to run the program, but
  ## absolutely need to change back to the original directory on
  ## exit.
  owd <- setwd(path)
  on.exit(setwd(owd))

  system(sprintf("%s -f .", program),
         ignore.stdout=!verbose, ignore.stderr=!verbose)
}

## Below here are internal functions

## Conversion between parameter formats for the reference
## implementation (faster-traitdiversity) and the tree data
## structures.

## Pending ugliness: What about the duplication of rho, s?

## Order and names of parameters in params.m and the core traits.  The
## order here is critical for both.  Note that "s" and "rho" are
## "seed_mass" and "wood_dens" (respectively) in params.m, but will
## have their names translated on read/write.
reference.params.names <- function()
  c("eta", "theta", "b", "a1", "B1", "a2", "B2", "a3",
    "a4", "B4", "n_area", "c_p1", "c_p2", "c_Rl", "c_Rs", "c_Rb",
    "c_Rr", "k_b", "k_r", "Y", "c_bio", "c_acc", "c_r1", "c_r2",
    "c_ext", "Pi_0", "c_d0", "c_d1", "c_d2", "c_d3", "s",
    "rho", "c_s0", "mean_disturbance_interval")
reference.traits.names <- function()
  c("lma", "hmat", "rho", "s")

## Output files we expect from running falster-growthmodel:
reference.output.files <- function()
  ## Files for which we will build a matrix or a data.frame.  The
  ## actual filename is <idx>_<type>.txt, with type from this vector.
  ## The first (popn) is special, and will be loaded as a data.frame.
  c("popn",
    "age", "coh_m", "coh_n", "coh_r", "bound_m", "bound_s", "bound_n",
    "bound_r", "d_bound_m", "d_bound_s", "d_bound_r")

## Generate a temporary transfer object (a small list) from a
## Parameters object, p.
reference.from.parameters <- function(p) {
  check.parameters.for.reference(p)

  ## Get the parameters from the first strategy:
  params <- p[[1]]$parameters[reference.params.names()]
  names(params) <- reference.params.names()

  ## These come from the overall simulation parameters:
  in.parameters <- c("c_ext", "Pi_0", "mean_disturbance_interval")
  params[in.parameters] <- p$parameters[in.parameters]

  ## Core traits:
  traits <- t(sapply(p$strategies, function(x)
                     unlist(x$parameters[reference.traits.names()])))

  ## Seed rain:
  seed.rain <- p$seed_rain

  ## And return a list
  list(params=params,
       traits=traits,
       seed.rain=seed.rain)
}

## Generate a Parameters object from the temporary transfer object (a
## small list).
parameters.from.reference <- function(obj) {
  check.transfer.object(obj)
  params <- obj$params
  traits <- obj$traits
  seed.rain <- obj$seed.rain

  in.parameters <- c("c_ext", "Pi_0", "mean_disturbance_interval")
  p <- new(Parameters)
  p$set_parameters(params[in.parameters])

  common <- params[setdiff(names(params), in.parameters)]
  for (i in seq_len(nrow(traits)))
    p$add_strategy(new(Strategy, modifyList(as.list(traits[i,]), common)))

  p$seed_rain <- seed.rain

  p
}

## Write the temporary transfer object out as the files "params.m" and
## "Stoch.txt" in some directory.
write.reference.parameters <- function(obj, path) {
  check.transfer.object(obj)
  dir.create(path, FALSE)

  ## 1: The params.m file:
  params <- obj$params
  ## Log the mean disturbance interval
  i <- match("mean_disturbance_interval", names(params))
  params[[i]] <- log10(params[[i]])
  names(params)[i] <- paste0("log_", names(params)[i])
  ## Rename variables as expected by falster-
  params <- rename(params, c(rho="wood_dens", s="seed_mass"))

  ## Lazy, but nice, double->string conversion at appropriate
  ## accuracy using %s.
  writeLines(c("function p=params",
               sprintf("p.%s=%s;", names(params), unlist(params))),
             file.path(path, "params.m"))

  ## 2: Stoch.txt
  traits <- obj$traits
  seed.rain <- obj$seed.rain
  n <- nrow(traits)
  tmp <- c(`%T`=0, NumRes=n, c(t(log10(traits))), seed.rain)
  names(tmp)[seq(3, length=n*4)] <-
    sprintf("y%d", seq(0, length=n*4))
  names(tmp)[seq(3+n*4, length.out=n)] <-
    sprintf("p%d", seq(0, length=n))
  write.table(data.frame(t(tmp), check.names=FALSE),
              file.path(path, "Stoch.txt"), sep="\t", row.names=FALSE,
              quote=FALSE)
}

## Read in the files "params.m" and "Stoch.txt" as a temporary
## transfer object.
read.reference.parameters <- function(path) {
  ## 1: Read the parameters
  params <- sub(";$", "", readLines(file.path(path, "params.m"))[-1])
  params <- do.call(rbind, strsplit(params, "="))
  params <- structure(as.list(as.numeric(params[,2])),
                      names=gsub("^p\\.", "", params[, 1]))
  ## Un-log the mean disturbance interval:
  i <- match("log_mean_disturbance_interval", names(params))
  params[[i]] <- 10^params[[i]]
  names(params)[i] <- sub("^log_", "", names(params)[i])
  ## Rename variables as expected by falster-traitdiversity:
  params <- rename(params, c(wood_dens="rho", seed_mass="s"))
  if (!identical(names(params), reference.params.names()))
    stop("Unexpected parameter names/order in params.m")

  ## 2: Core traits:
  core <- read.tsv.matrix(file.path(path, "Stoch.txt"), TRUE)

  ## Number of species/strategies:
  n <- core[[2]]

  ## The four core traits (translate *out* of log basis)
  traits <- matrix(10^core[seq(3, length.out=n * 4)], n, 4, byrow=TRUE)
  colnames(traits) <- c("lma", "hmat", "rho", "s")
  ## But ignore the specified values for wood density and seed mass,
  ## and use the version in params.m (see note within
  ## falster-traitdiversity in readme.md and in Strategy::set_traits).
  traits[, c("rho", "s")] <- rep(unlist(params[c("rho", "s")]), each=n)

  ## The seed rain:
  seed.rain <- core[seq(3 + n*4, length.out=n)]

  list(params=params,
       traits=traits,
       seed.rain=seed.rain)
}

## Utilities for checking that things are sensible:
check.transfer.object <- function(obj) {
  if (!identical(names(obj), c("params", "traits", "seed.rain")))
    stop("Invalid object")
  invisible(TRUE)
}

check.parameters.for.reference <- function(p) {
  if (p$size == 0) {
    stop("Cannot create parameters without strategy")
  } else if (p$size > 1) {
    ## Check that all strategies are the same aside from core traits:
    exclude <- function(x, ex) x[!(names(x) %in% ex)]
    pars <- sapply(p$strategies, function(x)
                   exclude(unlist(x$parameters), c("lma", "hmat")))
    err <- apply(pars, 1, diff) > 0
    if (any(err > 0))
      stop(sprintf("Parameters differ among stragies for %s",
                   paste(names(err[err > 0]), collapse=", ")))
  }
  invisible(TRUE)
}

check.reference.output <- function(n, path) {
  re <- "^([0-9]+)_age.txt$"
  path.output <- file.path(path, "res")
  strategy.idx <- as.integer(sub(re, "\\1", dir(path.output, pattern=re)))

  if (!isTRUE(all.equal(strategy.idx, seq(0, length=n))))
    stop("Incorrect number of strategies")

  idx <- seq(0, length.out=n)
  files.each <- reference.output.files()
  files <- sprintf("%d_%s.txt", idx, rep(files.each, each=n))
  other <- c("params.m", "patch_age.txt", "stand.txt", "Strategy.txt")
  expected <- c(files, other)

  ok <- file.exists(file.path(path.output, expected))
  if (any(!ok))
    stop(sprintf("Could not find file(s) in %s: %s", path.output,
                 paste(expected[!ok], collapse=", ")))

  extra <- setdiff(dir(path.output), expected)
  if (length(extra) > 0)
    stop(sprintf("Additional files found in %s: %s", path.output,
                 paste(extra, collapse=", ")))
  invisible(TRUE)
}

## General utility -- might be useful elsewhere.
rename <- function(obj, tr) {
  i <- match(names(tr), names(obj))
  if (any(is.na(i)))
    stop("Expected names missing from obj: ",
         paste(names(tr)[is.na(i)], collapse=", "))
  names(obj)[i] <- tr
  obj
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
