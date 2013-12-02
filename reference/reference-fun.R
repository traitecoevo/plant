## These utilities are here only temporarily, to streamline the
## 'reference.R' file as much as possible.  Once everything is
## working, these will probably get incorporated into the reference
## loading functions.

## This should perhaps be run in the reference comparison functions,
## then I'd never have to bother with doing the translation for the
## actual model output.  Plus I will have the strategies on hand when
## doing this.
tidyup.reference <- function(output, strategy) {
  values <- list(
    mass.leaf         = output$strategies[[1]]$bound_m,
    survival          = output$strategies[[1]]$bound_s,
    seeds             = output$strategies[[1]]$bound_r,
    density.mass.leaf = output$strategies[[1]]$bound_n)
  values$height <-
    translate(values$mass.leaf, mass.leaf.to.height, strategy)
  values$density <-
    values$mass.leaf * dmldh(values$height, strategy)
  values$log.density <- log(values$density)
  values$mortality <- -log(values$survival)
  ## Then order these so that the first four are the same order as
  ## CohortTop produces, and the rest are hopefully sensible.
  values <- values[c("height", "mortality", "seeds", "log.density",
                     "mass.leaf", "survival",
                     "density", "density.mass.leaf")]

  ## We don't have access to the rate in change in density(mass.leaf)
  ## over time; not sure why, but it looks like that is never reported
  ## from the EBT.
  rates <- list(
    mass.leaf         = output$strategies[[1]]$d_bound_m,
    survival          = output$strategies[[1]]$d_bound_s,
    seeds             = output$strategies[[1]]$d_bound_r,
    log.density       = NULL)

  ## Not checked:
  ## d(survival)/dt -> exp(-(d(mortality)/dt / mortality))
  ##   -- D[-Log[m[t]], t]
  ## d(height)/dt -> dm/dt * dh/dm
  ##   -- D[a (m[t]/lma) ^ B1, t]
  rates$mortality <-
    -1 / values$survival * rates$survival
  rates$height <-
    rates$mass.leaf * dhdml(values$mass.leaf, strategy)

  i.schedule <- grep("^add.cohort.[0-9]+$", names(output$patch_age))
  schedule <- as.matrix(output$patch_age[i.schedule])
  dimnames(schedule) <- NULL

  list(time=output$patch_age$age,
       values=values,
       rates=rates,
       light.env=output$light.env,
       schedule=schedule)
}

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

## This is more generally useful?
translate <- function(m, f, ...)
  array(f(m, ...), dim(m))

## This has been copied around enough to be useful.  Perhaps change
## this to do 'is.valid.for.reference' followed by a more general
## compare function?
compare.parameters <- function(p1, p2)
  isTRUE(all.equal(tree:::reference.from.parameters(p1),
                   tree:::reference.from.parameters(p2)))

## Collect all variables at each cohort introduction.
run.ebt <- function(ebt, max.t=Inf) {
  pad.matrix <- function(x) {
    n <- max(sapply(x, length))
    t(sapply(x, function(i) c(i, rep(NA, n-length(i)))))
  }
  collect <- function() {
    values <- matrix(ebt$ode_values, 4)
    rates  <- matrix(ebt$ode_rates,  4)
    rownames(values) <- rownames(rates) <-
      c("height", "mortality", "seeds", "log.density")
    light.env <- ebt$patch$environment$light_environment$xy
    colnames(light.env) <- c("height", "canopy.openness")
    list(time=ebt$time,
         values=values,
         rates=rates,
         light.env=light.env)
  }
  combine <- function(res) {
    time <- sapply(res, "[[", "time")

    values <- lapply(res, "[[", "values")
    rates  <- lapply(res, "[[", "rates")
    vars <- rownames(values[[1]])
    values <- lapply(vars, function(v)
                     pad.matrix(lapply(values, function(x)
                                       unname(x[v,]))))
    rates <- lapply(vars, function(v)
                    pad.matrix(lapply(rates, function(x)
                                      unname(x[v,]))))
    names(values) <- names(rates) <- vars

    light.env <- lapply(res, "[[", "light.env")

    list(time=time,
         values=values,
         rates=rates,
         light.env=light.env)
  }

  res <- list()
  ebt$reset()

  while (ebt$cohort_schedule$remaining > 0 && ebt$time < max.t) {
    ebt$run_next()
    res <- c(res, list(collect()))
  }

  combine(res)
}

## Resample the reference output down to a given set of times:
resample.reference <- function(ref) {
  idx <- c(which(apply(ref$schedule, 1, any))[-1], nrow(ref$schedule))
  sample <- function(y.old)
    y.old[idx,,drop=FALSE]
  ret <- list(time  = ref$time[idx],
              values= lapply(ref$values, sample),
              rates = lapply(ref$rates,  sample),
              light.env = ref$light.env[idx])
  remove.boundary.cohort(ret)
}

rescale.reference <- function(ref, t.new) {
  t.old <- ref$time
  scale <- function(y.old) {
    nok <- length(na.omit(y.old))
    t.min <- t.old[!is.na(y.old)][1]
    if (nok < 2)
      ret <- numeric(0)
    else if (length(na.omit(y.old)) < 3)
      ret <- approx(t.old, y.old, xout=t.new[t.new >= t.min])$y
    else
      ret <- spline(t.old, y.old, xout=t.new[t.new >= t.min])$y
    c(rep(NA, length(t.new) - length(ret)), ret)
  }
  f <- function(m)
    apply(m, 2, scale)
  ret <- list(time  = t.new,
              values= lapply(ref$values, f),
              rates = lapply(ref$rates,  f))
  remove.boundary.cohort(ret)
}

remove.boundary.cohort <- function(ref) {
  ## Do we always drop the last two columns?  It might depend on if a
  ## cohort is introduced.
  f <- function(m) {
    if (is.null(m))
      return(m)
    n <- nrow(m)
    if (ncol(m) != n + 2)
      stop("Unexpected size for output matrix")
    m <- m[,-(n + 1:2)]
    m[cbind(1:(n-1), 2:n)] <- NA
    m
  }

  ref$values <- lapply(ref$values, f)
  ref$rates  <- lapply(ref$rates,  f)
  ref
}

## We expect to see 'evolve' in the directory pointed at by
## PATH_EVOLVE:
##   1. exists
##   2. contains the 'evolve' program in src
##   3. has the correct git hash
##
## This is a bit of a fragile hack around the fact that nobody will
## actually have this installed globally on their system (in PATH or
## something).  Once falster-traitdiversity is opened up on github, we
## could just download it into tree's installed directory as needed.
locate.evolve <- function() {
  path.evolve <- Sys.getenv("PATH_EVOLVE")
  path.evolve.cmd <- file.path(path.evolve, "src")
  if (!(file.exists(path.evolve) &&
        file.exists(file.path(path.evolve.cmd, "evolve"))))
    stop("Could not find 'evolve': expected in ", path.evolve.cmd)

  ## Executed in a subshell, so no risk of changing directory.
  evolve.version <-
    system(sprintf("(cd %s && git rev-parse HEAD)", path.evolve),
           intern=TRUE)
  if (evolve.version != "6e63e505f37827303346f44c5886715bd080fd2c")
    warning("'evolve' version has changed")
  path.evolve.cmd
}
