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
##' \item{patch_age}{\code{data.frame} with columns \code{age} (age in
##' years), \code{density} (density of patches of that age),
##' \code{gamma} (parameter of the disturbance function).  After that
##' are columns \code{add.cohort.1}, \code{add.cohort.2}, etc, which
##' are logical vectors indicating if a cohort for a strategy was
##' added at that time point.}
##' \item{seed_rain_out}{Vector of output seed rain, one entry per
##' strategy/species}
##' \item{light_env}{Light environment as one matrix per time point.
##' Each matrix has columns \code{height} (height) and
##' \code{canopy.openness}.  These are always 50 rows long.}
##' \item{popn}{A list of lists, one per strategy.  Each strategy's
##' list has components:
##' \describe{
##'   \item{age}{Age of the cohort at a given point in time.}
##'   \item{coh_m}{Leaf mass of the mean individual of a cohort}
##'   \item{coh_n}{Mean density (leaf mass, see eqn 22)}
##'   \item{coh_r}{Mean seed production}
##'   \item{bound_m}{Leaf mass of the boundary individual of a cohort}
##'   \item{bound_s}{Boundary seed production}
##'   \item{bound_n}{Boundary density (leaf mass, see eqn 22)}
##'   \item{bound_r}{Boundary seed production}
##'   \item{d_bound_m}{Time derivative of \code{bound_m}}
##'   \item{d_bound_s}{Time derivative of \code{bound_s}}
##'   \item{d_bound_r}{time derivative of \code{bound_r}}
##'   \item{popn}{Matrix of physiological/demographic summary
##'      statistics, one row per time.}
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
  traits.output <- popn$strategy[, sprintf("Tr%d", 0:3), drop=FALSE]
  ok <- isTRUE(all.equal(unname(traits.input), unname(traits.output),
                         tolerance=1e-5))
  if (!ok)
    stop("Input and output traits appear to vary")
  if (!isTRUE(all.equal(unname(popn$strategy[, "X_start"]),
                        unname(p$seed_rain), tolerance=1e-5)))
    stop("Input and output seed_rain appear to vary")

  ## Extract the output seed rain from strategy, and drop the rest of
  ## the strategy element, as it is just input parameters.
  popn$seed_rain_out <- popn$strategy[, "X_end"]
  popn <- popn[names(popn) != "strategy"]

  ## Process the patch_age element to separate out the
  ## age/density/disturbance part (as a data.frame) and the light
  ## environment part (as a list of matrices)
  tmp <- popn$patch_age
  i.h <- grep("^H[0-9]+$", colnames(tmp))
  i.e <- grep("^E[0-9]+$", colnames(tmp))
  light.env <- lapply(seq_len(nrow(tmp)), function(i)
                      cbind(height=unname(tmp[i, i.h]),
                            canopy.openness=unname(tmp[i, i.e])))

  ## Alternatively, but possibly more confusingly:
  ## n.steps <- nrow(tmp)
  ## n.heights <- length(i.e)
  ## array(tmp[,c(i.h, i.e)], c(n.steps, n.heights, 2),
  ##              dimnames=list(step=NULL, height=NULL,
  ##                c("height", "canopy.openness")))

  popn$patch_age <- data.frame(age=tmp[, "Age"],
                               density=tmp[, "p(a,t)"],
                               gamma=tmp[, "gamma(a)"])

  ## Load a single strategy object, with index 'idx':
  load.strategy <- function(idx) {
    fmt <- sprintf("%d_%%s.txt", idx)
    ret <- lapply(reference.output.files(), function(x)
                  read.tsv.matrix(file.path(path.output, sprintf(fmt, x)),
                                  x == "popn"))
    names(ret) <- reference.output.files()
    for (i in seq_along(ret)[-1]) { # skip $popn
      x <- ret[[i]]
      if (!all(is.na(x[,ncol(x)])))
        stop("Unexpected non-NA input in final column")
      ret[[i]] <- x[,-ncol(x)]
    }
    ret
  }
  popn$strategies <- lapply(seq(0, length.out=p$size), load.strategy)
  popn$light.env  <- light.env
  popn$parameters <- p

  ## Add the cohort introduction information to "age":
  ##
  ## NOTE: x$popn[,"NoCohorts"] turns out not to be correct with
  ## respect to the boundaries.  It is correct for the mean cohorts
  ## though (so I guess this is basically a fencepost error).
  ## Manually pull out counts for both this way:
  introduction <- function(y)
    diff(c(0, rowSums(!is.na(y)))) > 0
  popn$introduction.bound <- sapply(popn$strategies, function(x)
                                    introduction(x$bound_m))
  popn$introduction.mean <- sapply(popn$strategies, function(x)
                                   introduction(x$coh_m))

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
##' @param verbose Logical indicating if the model output should be
##' printed (TRUE by default).
##' @param evolve Path to the "evolve" program.  If not specified,
##' "evolve" is assumed to be installed within tree by
##' \code{install.evolve}
##' @author Rich FitzJohn
##' @rdname run.reference
##' @export
run.reference <- function(path, p=NULL, verbose=TRUE, evolve=NULL) {
  if (is.null(evolve)) {
    if (!evolve.is.installed())
      stop("'evolve' not installed - run install.evolve() first")
    evolve <- file.path(path.evolve.dir(), "evolve")
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

  system(sprintf("%s -f .", evolve),
         ignore.stdout=!verbose, ignore.stderr=!verbose)
}

##' @export
##' @rdname run.reference
install.evolve <- function(reinstall=FALSE, verbose=FALSE,
                           url=NULL) {
  path <- path.evolve.dir()
  if (file.exists(path)) {
    if (reinstall)
      uninstall.evolve()
    else
      return(invisible(FALSE))
  }

  if (is.null(url)) {
    if (file.exists("~/.falster-traitdiversity_path"))
      url <- readLines("~/.falster-traitdiversity_path")
    else
      url <- "git@github.com:dfalster/Falster-TraitDiversity.git"
  }
  hash <- "6e63e505f37827303346f44c5886715bd080fd2c"
  intern <- !verbose

  system(paste("git clone", url, path), intern=intern)
  system(sprintf(paste("git --git-dir=%s/.git --work-tree=%s",
                       "checkout -b tree_version %s"),
                 path, path, hash), intern=intern)
  system(sprintf("cd %s/src && make", path), intern=intern)
  file.copy(file.path(path, "src", "evolve"), path)

  invisible(TRUE)
}

##' @export
##' @rdname run.reference
uninstall.evolve <- function() {
  unlink(path.evolve.dir(), recursive=TRUE)
}

##' @export
##' @rdname run.reference
evolve.is.installed <- function()
  file.exists(file.path(path.evolve.dir(), "evolve"))

path.evolve.dir <- function()
  file.path(path.package("tree"), "falster-traitdiversity")

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
  in.parameters <- c("c_ext", "Pi_0")
  params[in.parameters] <- p$parameters[in.parameters]
  params$mean_disturbance_interval <- p$disturbance$mean_interval

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

reference.compare.parameters <- function(p1, p2)
  isTRUE(all.equal(reference.from.parameters(p1),
                   reference.from.parameters(p2)))

## Generate a Parameters object from the temporary transfer object (a
## small list).
parameters.from.reference <- function(obj) {
  check.transfer.object(obj)
  params <- obj$params
  traits <- obj$traits
  seed.rain <- obj$seed.rain

  in.parameters <- c("c_ext", "Pi_0")
  p <- new(Parameters)
  p$set_parameters(params[in.parameters])
  p$disturbance <- new(Disturbance, params[["mean_disturbance_interval"]])

  common <- params[setdiff(names(params),
                           c(in.parameters, "mean_disturbance_interval"))]
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
  tmp <- c(T=0, NumRes=n, c(t(log10(traits))), seed.rain)
  names(tmp)[seq(3, length=n*4)] <-
    sprintf("y%d", seq(0, length=n*4))
  names(tmp)[seq(3+n*4, length.out=n)] <-
    sprintf("p%d", seq(0, length=n))
  write.tsv.matrix(data.frame(t(tmp), check.names=FALSE),
                   file.path(path, "Stoch.txt"))
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
  if (with.header) {
    h <- scan(con, nlines=1, quiet=TRUE, what=character())
    h[[1]] <- sub("^%", "", h[[1]])
  }

  d1 <- scan(con, nlines=1, quiet=TRUE, sep="")
  d2 <- scan(con, quiet=TRUE, sep="")
  m <- matrix(c(d1, d2), ncol=length(d1), byrow=TRUE)

  if (with.header)
    colnames(m) <- h

  m[is.nan(m)] <- NA
  m
}

write.tsv.matrix <- function(x, file) {
  if (substr(names(x)[[1]], 1, 1) != "%")
    names(x)[[1]] <- paste0("%", names(x)[[1]])
  write.table(x, file, sep="\t", row.names=FALSE, quote=FALSE)
}

##' @export
##' @param output Reference output from \code{load.reference.output}
##' @rdname load.reference.parameters
tidyup.reference.output <- function(output) {
  ## First element is time - that's easy to get.
  time <- output$patch_age$age
  ## The light environment is also fairly easy to get.
  light.env <- output$light.env

  output$patch_age$pr.surv <-
    output$patch_age$density / output$patch_age$density[[1]]

  ## The things we *really* want are height, log.mortality, seeds,
  ## log.density (*height*) and pr.survival.birth, to match the output
  ## of tree.
  f <- function(i) {
    x <- output$strategies[[i]]
    p <- output$parameters
    strategy <- p[[i]]

    values <- list(mass.leaf         = x$bound_m,
                   survival          = x$bound_s,
                   seeds             = x$bound_r,
                   density.mass.leaf = x$bound_n)

    ## Extracting the per-cohort birth probability of survival is
    ## tricky.
    pr.surv.birth <- output$patch_age$pr.surv[output$introduction.bound[,i]]
    if (length(pr.surv.birth) != ncol(values$mass.leaf))
      stop("Unexpected dimension sizes")
    tmp <- array(rep(pr.surv.birth, each=nrow(values$mass.leaf)),
                 dim(values$mass.leaf))
    tmp[is.na(values$mass.leaf)] <- NA
    values$pr.survival.birth <- tmp

    values$log.mortality <- -log(values$survival)
    values$height <- mass.leaf.to.height(values$mass.leaf, strategy)
    ## Log density of *height* - see EBT.md, eq:n_conversion2.
    ##
    ## I have an awful feeling that we set the initial conditions
    ## differently between implementations though.  In tree, it is set
    ## to
    ##   seed_rain * pr_germ / (dh/dt)
    ## but I'm not sure what it's set to in faster-traitdiversity.
    values$log.density <- log(values$density.mass.leaf *
                              dmldh(values$height, strategy))

    values$seeds <- values$seeds / p$parameters$Pi_0

    v <- c("height", "log.mortality", "seeds", "log.density",
           "pr.survival.birth")
    ret <- aperm(list.to.array(values[v]), c(3, 1, 2))
  }

  species <- lapply(seq_along(output$strategies), f)

  list(time=time, species=species, light.env=light.env,
       schedule=output$introduction.bound)
}

##' @export
##' @param output Reference output from \code{tidyup.reference.output}
##' @rdname load.reference.parameters
resample.reference.output <- function(output) {
  idx <- c(which(apply(output$schedule, 1, any)), nrow(output$schedule))
  ## Subset:
  output$time      <- output$time[idx]
  output$species   <- lapply(output$species, function(x) x[,idx,,drop=FALSE])
  output$light.env <- output$light.env[idx]
  output$schedule  <- NULL
  output
}

##' @export
##' @param output Reference output from \code{tidyup.reference.output}
##' @rdname load.reference.parameters
schedule.from.reference <- function(output) {
  n <- length(output$species)
  sched <- new(CohortSchedule, n)
  for (i in seq_len(n))
    sched$set_times(output$time[output$schedule[,i]], i)
  sched$max_time <- max(output$time)
  sched$ode_times <- output$time
  sched
}

## Utilities used -- these need documenting and testing, probably.
## Conversion functions.  Where should these go so that they are
## fairly efficient to use?  Going through Plant is slow for these.
mass.leaf.to.height <- function(mass.leaf, strategy) {
  pars <- strategy$parameters
  a1 <- pars$a1
  B1 <- pars$B1
  lma <- pars$lma
  a1 * (mass.leaf / lma) ^ B1
}
height.to.mass.leaf <- function(height, strategy) {
  pars <- strategy$parameters
  a1 <- pars$a1
  B1 <- pars$B1
  lma <- pars$lma
  (height / a1) ^ (1 / B1) * lma
}
dhdml <- function(mass.leaf, strategy) {
  pars <- strategy$parameters
  a1 <- pars$a1
  B1 <- pars$B1
  lma <- pars$lma
  a1 * B1 / lma * (mass.leaf / lma) ^ (B1 - 1)
}
dmldh <- function(h, strategy) {
  pars <- strategy$parameters
  a1 <- pars$a1
  B1 <- pars$B1
  lma <- pars$lma
  lma / (a1 * B1) * (h / a1) ^ (1 / B1 - 1)
}

## Disturbance cumulative probability of survival; used for working
## out how long to run the simulation for.  This is of unknown origin
## within the falster-traitdiversity, but is the value that
## corresponds to
##   2.633*a_mean /3.0*4.0;
## within src/base/ebt/site.cpp, site::solve_patchage_dist()
##
## It's here so that we can do
##   disturbance$cdf(tree:::refrence.pr.survival.eps)
## and get a running time that agrees with falster-traitdiversity.
reference.pr.survival.eps <- 6.25302620663814e-05
