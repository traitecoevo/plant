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

##' Sequence in log space
##'
##' Unlike the billions of options for \code{seq}, only
##' \code{length.out} is supported here, and both \code{from} and
##' \code{to} must be provided.
##' @title Sequence in Log Space
##' @param from Starting point
##' @param to Ending point
##' @param length.out Number of points to generate
##' @author Rich FitzJohn
##' @export
seq_log <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out=length.out))
}

##' @export
##' @param r range (i.e., c(from, to)
##' @rdname seq_log
seq_log_range <- function(r, n) {
  seq_log(r[[1]], r[[2]], n)
}


collect <- function(k, x, empty=list(), each=identity, after=identity,
                    loop=sapply) {
  if (length(x) == 0) {
    empty
  } else {
    after(loop(x, function(el) each(el[[k]])))
  }
}

rbind_list <- function(x) {
  do.call(rbind, as.list(x))
}

#' Calculate the gradient of a function by numerical approximation
#'
#' The function ‘gradient_fd’ calculates a numerical approximation
#' of the first derivative of ‘func’ at the point ‘x’ using finite
#' differencing. It is assumed ‘func’ is a scalar value function.
#' @param func a function with a scalar real result
#' @param  x a real scalar or vector argument to func, indicating
#' the point(s) at which the gradient is to be calculated.
#' @param  dx Interval over which derivative is calculated
#' @param  log_scale=TRUE Determines whether derivative is taken
#' with respect to raw or log-transformed x values. The latter is
#' equiavlent to taking the derivative of a function
#' \code{g(f(exp(x)))} with respect to log x and is useful when
#' x is log-normally distributed.
#' @author Daniel Falster
#' @export
gradient_fd <- function(func, x, dx, log_scale=FALSE) {
  if (!log_scale) {
    xx <- x + c(-1, 1) * dx/2
  } else {
    xx <- x * exp(c(-1, 1) * dx/2)
  }
  yy <- func(xx)
  (yy[[2]] - yy[[1]]) / dx
}

##' @export
git_sha <- function() {
  system("git rev-parse --short HEAD", intern=TRUE)
}

splinefun_log <- function(x, y, ...) {
  f <- splinefun(log(x), y, ...)
  function(x) {
    f(log(x))
  }
}

## Really simple rejection sampling, assuming:
##   1. univariate function
##   2. uniform approximation to f is a reasonable upper bound
##   3. that f_max > f for all x
rejection_sample <- function(n, f, bounds, f_max=NULL, log_space=TRUE,
                             action_nopositive=message) {
  rejection_sample_iter <- function() {
    x <- runif(n, bounds[[1]], bounds[[2]])
    fx <- f(x)
    keep <- fx / f_max
    x[runif(n) < keep]
  }

  if (log_space) {
    bounds <- log(bounds)
    f_orig <- f
    f <- function(x) f_orig(exp(x))
  }

  ## Hard coded, and possibly not very clever:
  if (is.null(f_max)) {
    f_max <- max(f(seq(bounds[[1]], bounds[[2]], length.out=501))) * 1.2
  }
  res <- numeric(0)
  if (f_max <= 0) {
    action_nopositive("No positive values found")
  } else {
    while (length(res) < n) {
      res <- c(res, rejection_sample_iter())
    }
    res <- res[seq_len(n)]
    if (log_space) {
      res <- exp(res)
    }
  }
  res
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

has_attr <- function(x, which, exact=TRUE) {
  !is.null(attr(x, which, exact))
}
