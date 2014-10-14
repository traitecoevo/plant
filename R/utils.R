list_to_array <- function(x) {
  if (length(unique(lapply(x, dim))) > 1)
    stop("More than one dimension")

  dimnames2 <- function(x) {
    dn <- dimnames(x)
    if (is.null(dn)) rep(list(NULL), length(dim(x))) else dn
  }
  array(unlist(x),
        c(dim(x[[1]]), length(x)),
        dimnames=c(dimnames2(x[[1]]), list(names(x))))
}

pad_matrix <- function(x) {
  if (is.matrix(x[[1]])) {
    nc <- max(sapply(x, ncol))
    nr <- nrow(x[[1]])
    lapply(x, function(i) cbind(i, matrix(NA, nr, nc - ncol(i))))
  } else {
    nc <- max(sapply(x, length))
    lapply(x, function(i) c(i, rep(NA, nc - length(i))))
  }
}

pad_list_to_array <- function(x) {
  list_to_array(pad_matrix(x))
}

##' Get last element from an object
##' @title Get Last Element
##' @param x An object that can be subset with \code{[[}
##' @author Rich FitzJohn
##' @export
last <- function(x)
  x[[length(x)]]
##' @rdname last
##' @export
`last<-` <- function(x, value) {
  x[[length(x)]] <- value
  x
}

##' ##' @rdname last
##' @export
first <- function(x)
  x[[1]]
##' @rdname last
##' @export
`first<-` <- function(x, value) {
  x[[1]] <- value
  x
}

##' Runs an expression and caches the result in a file.  Very basic,
##' no checks.
##'
##' @title Run and Cache an Expression
##' @param expr An expression to be evaluated
##' @param filename A filename to store the result (will be saved via
##' \code{\link{saveRDS}}
##' @param regenerate Logical value indicating if results should be
##' regenerated.
##' @return The result of evaluating \code{expr}, possibly loaded from
##' \code{filename} rather than being rerun.
##' @author Rich FitzJohn
##' @export
run_cached <- function(expr, filename, regenerate=FALSE) {
  if (file.exists(filename) && !regenerate) {
    res <- readRDS(filename)
  } else {
    res <- eval.parent(substitute(expr))
    saveRDS(res, filename)
  }
  res
}

##' Add points and lines for a cobweb plot showing convergence to an
##' equilibrium.
##'
##' Assumed (unchecked) that there are two columns corresponding to
##' "from" and "to.
##' @title Cobweb Plot
##' @param m Two column matrix
##' @param ... Additional parameters passed to \code{lines}
##' @author Rich FitzJohn
##' @export
cobweb <- function(m, ...) {
  lines(rep(m[,1], each=2), c(t(m)), ...)
}

git_sha <- function(package=.packageName) {
  read_system_file_if_exists <- function(path) {
    filename <- system.file(path, package=package)
    if (filename == "") character(0) else readLines(filename)
  }
  sha    <- read_system_file_if_exists("git/sha")
  status <- read_system_file_if_exists("git/status")
  if (length(status) > 0) {
    attr(sha, "status") <- parse_git_status(status)
  }
  sha
}

parse_git_status <- function(x) {
  cbind(index=substr(x, 1, 1),
        work=substr(x, 2, 2),
        path=substr(x, 4, nchar(x)))
}


## This does all sorts of nasty, but will be useful for development
## from R without the hellish reload cycle of Rcpp modules.  Using
## ':::' will then *not* work.
reload_r <- function(path=.TREE_PATH) {
  files <- dir(file.path(path, "R"), pattern=glob2rx("*.R"),
               full=TRUE)
  for (f in files) {
    source(f, local=FALSE)
  }
  invisible(TRUE)
}

recycle_simple <- function(x, n) {
  if (length(x) == 1) {
    x <- rep(x, length=n)
  } else if (length(x) != n) {
    stop("Length must be 1 or ", n)
  }
  x
}

##' Thin wrapper around \code{nleqslv} and \code{dfsane}
##' @title Thin wrapper around nleqslv and dfsane
##' @param x Starting point
##' @param fn Function to solve
##' @param tol Tolerance (for \code{nleqslv} this will be both
##' absolute and relative)
##' @param maxit Maximum number of iterations.  The number of function
##' evaluations will likely exceed this.
##' @param solver The solver to use.  Either "nleqslv" or "dfsane" for now.
##' @export
nlsolve <- function(x, fn, tol=1e-6, maxit=100, solver="nleqslv") {
  solver <- match.arg(solver, c("nleqslv", "dfsane"))

  if (solver == "nleqslv") {
    res <- nlsolve_nleqslv(x, fn, tol, maxit)
  } else if (solver == "dfsane") {
    res <- nlsolve_dfsane(x, fn, tol, maxit)
  } else {
    stop("Unknown solver ", solver)
  }

  if (attr(res, "failed")) {
    stop(sprintf("Solver has likely failed: code=%d, msg: %s",
                 attr(res, "code"), attr(res, "message")),
         immediate.=TRUE)
  }

  res
}

nlsolve_nleqslv <- function(x, fn, tol=1e-6, maxit=100) {
  control <- list(xtol=tol, ftol=tol, maxit=maxit)
  sol <- nleqslv::nleqslv(x, fn, global="none", control=control)
  code <- sol$termcd
  res <- sol$x
  attributes(res) <- nlsolve_nleqslv_attr(sol)
  res
}

nlsolve_nleqslv_attr <- function(sol) {
  list(y=sol$fvec, # different to dfsane
       iter=sol$iter,
       feval=sol$nfcnt, # does not include jacobian evals
       code=sol$termcd,
       message=sol$message,
       failed=sol$termcd > 2 || sol$termcd < 0,
       solver="nleqslv")
}

nlsolve_dfsane <- function(x, fn, tol=1e-6, maxit=100) {
  control <- list(tol=tol, maxit=maxit, trace=FALSE)
  sol <- BB::dfsane(x, fn, control=control, quiet=TRUE)
  res <- sol$par
  attributes(res) <- nlsolve_dfsane_attr(sol)
  res
}

nlsolve_dfsane_attr <- function(sol) {
  list(y=sol$residual,
       iter=sol$iter,
       feval=sol$feval,
       code=sol$convergence,
       message=sol$message,
       failed=sol$convergence != 0,
       solver="dfsane")
}

failed <- function(x) {
  inherits(x, "try-error")
}

## Parallel version of Richardson extrapolation that will take
## advantage of tree's preference for computing multiple mutants at
## once.  Not yet set up to also do multiple x points at once.
##
## Defaults are set to match numDeriv::grad.
##' @export
gradient <- function(f, x,
                     eps=1e-4, d=0.0001, r=4, v=2,
                     zero_tol=sqrt(.Machine$double.eps/7e-7)) {
  ## This is not actually hard to relax, I just don't need that right
  ## now.  The big things are:
  ##   - recycle d appopriately
  ##   - x_up, x_down are matrices and are rbound/cbound together and
  ##     we require that from f
  ##   - once y is computed everything is the same.
  if (length(x) != 1) {
    stop("Requires scalar x")
  }
  ## Initial offset:
  h0 <- abs(d * x) + eps * (abs(x) < zero_tol)
  h <- h0 / (v^(seq_len(r) - 1))
  x_up   <- x + h
  x_down <- x - h

  ## This is the bit I'm trying to get: *all* the x points are
  ## computed at once.
  y <- matrix(f(c(x_up, x_down)), r, 2)

  ## Then bookkeeping galore.
  a <- (y[,1] - y[,2]) / (2 * h)
  m <- 1L
  while (m < r) {
    four_m <- 4 ^ m
    a_next <- numeric(r - m)
    for (i in seq_along(a_next)) {
      a_next[i] <- (a[i + 1L] * four_m - a[i])/(four_m - 1.0)
    }
    a <- a_next
    m <- m + 1L
  }
  a
}
