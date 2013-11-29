build.schedule <- function(p, nsteps, n.t, t.max, eps,
                           progress=FALSE, verbose=FALSE) {
  ebt <- new(EBT, p)
  ebt$cohort_schedule <- default.schedule(n.t, t.max)
  times <- ebt$cohort_schedule$times(1)

  history <- list()
  for (i in seq_len(nsteps)) {
    err <- run.with.lai(times, ebt)
    split <- err$total > eps
    times.new <- split.times(times, split)
    if (progress)
      history <- c(history, list(list(times=times, split=which(split),
                                      err=err, w=ebt$fitness(1))))
    times <- times.new
    if (!any(split, na.rm=TRUE))
      break
    if (verbose)
      message(sprintf("%d: Splitting %d times (%d)",
                      i, sum(split), length(times)))
  }
  if (progress)
    attr(times, "progress") <- history
  times
}

## This gets the leaf area error out at each cohort introduction;
## that's the error criterion we care about.
run.with.lai <- function(times, ebt) {
  idx <- 1
  ebt$reset()
  ebt$set_times(times, idx)
  err.lai <- list()
  while (ebt$cohort_schedule$remaining > 0) {
    ebt$run_next()
    err.lai <- c(err.lai, list(ebt$leaf_area_error(idx)))
  }
  err.lai <- do.call(rbind, pad.matrix(err.lai))
  err.w <- ebt$fitness_error(idx)
  total <- suppressWarnings(apply(rbind(err.lai, err.w), 2, max, na.rm=TRUE))

  list(lai=err.lai, w=err.w, total=total)
}

split.times <- function(times, i) {
  ## Upwind splitting scheme only, which means that we will never
  ## split the last interval [assuming OK for now].  Inefficiently
  ## interleaved with sort().  These issues can change easily enough
  ## later.  The aim is making sure that we don't introduce the same
  ## point twice; one from upstream and one from downstream.
  dt <- diff(times)
  i <- which(i)
  sort(c(times, times[i] - dt[i-1]/2))
}

## Will probably update these so that n.t *or* times can be given
build.schedule.fitness <- function(p, nsteps, n.t, t.max, cores=1,
                                   progress=FALSE, verbose=FALSE) {
  ebt <- new(EBT, p)
  ebt$cohort_schedule <- default.schedule(n.t, t.max)
  times <- ebt$cohort_schedule$times(1)

  w <- run.with.times(times, ebt)
  obj <- list(times=times, w=w)
  history <- list(list(obj))

  for (i in seq_len(nsteps)) {
    obj <- build.schedule.refine(obj$times, obj$w, ebt,
                                 cores=cores, verbose=verbose)
    history <- c(history, list(obj))
  }

  times <- obj$times
  if (progress)
    attr(times, "progress") <- history
  times
}

build.schedule.refine <- function(times, w, ebt, cores=1,
                                  verbose=FALSE) {
  run.with.insert <- function(i)
    run.with.times(insert.time(i, times), ebt)
  insert.time <- function(i, t) {
    j <- seq_len(i)
    c(t[j], (t[i] + t[i+1])/2, t[-j])
  }

  res <- unlist(mclapply(seq_len(length(times) - 1), run.with.insert,
                         mc.cores=cores))
  j <- which.max(abs(res - w))
  if (verbose)
    message(sprintf("Inserting time at %2.5f [%d]; fitness = %2.5f",
                    mean(times[j:(j+1)]), j, res[j]))

  list(times=insert.time(j, times), idx=j+1, w=res[j])
}

run.with.times <- function(times, ebt) {
  idx <- 1
  ebt$reset()
  ebt$set_times(times, idx)
  ebt$run()
  ebt$fitness(idx)
}

## Generally useful function...
last <- function(x)
  x[[length(x)]]

## Ocassionally useful...
pad.matrix <- function(x) {
  if (is.matrix(x[[1]])) {
    nc <- max(sapply(x, ncol))
    nr <- nrow(x[[1]])
    lapply(x, function(i) cbind(i, matrix(NA, nr, nc - ncol(i))))
  } else {
    nc <- max(sapply(x, length))
    lapply(x, function(i) c(i, rep(NA, nc - length(i))))
  }
}

run.cached <- function(expr, filename, regenerate=FALSE) {
  if (file.exists(filename) && !regenerate) {
    res <- readRDS(filename)
  } else {
    res <- eval.parent(substitute(expr))
    saveRDS(res, filename)
  }
  res
}
