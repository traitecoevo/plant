first <- function(x) {
  x[[1]]
}
`first<-` <- function(x, value) {
  x[[1]] <- value
  x
}

second <- function(x) {
  x[[2]]
}

last <- function(x) {
  x[[length(x)]]
}
`last<-` <- function(x, value) {
  x[[length(x)]] <- value
  x
}

list_to_array <- function(x) {
  if (length(unique(lapply(x, dim))) > 1L) {
    stop("More than one dimension")
  }

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

##' Sequence in log space
##'
##' Unlike the billions of options for \code{seq}, only
##' \code{length.out} is supported here, and both \code{from} and
##' \code{to} must be provided.  For completeness, \code{seq_range}
##' generates a range in non-log space.
##' @title Sequence in log space
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
seq_log_range <- function(r, length.out) {
  seq_log(r[[1]], r[[2]], length.out)
}

##' @export
##' @rdname seq_log
seq_range <- function(r, length.out) {
  seq(r[[1]], r[[2]], length.out=length.out)
}

validate <- function(x, ...) {
  UseMethod("validate")
}
##' @export
validate.FFW16_Parameters <- function(x, ...) {
  FFW16_Parameters__vdor(x)
}

loop <- function(X, FUN, ..., parallel=FALSE) {
  if (parallel) {
    parallel::mclapply(X, FUN, ..., mc.preschedule=FALSE)
  } else {
    lapply(X, FUN, ...)
  }
}

make_message_verbose <- function(verbose) {
  if (verbose) {
    function(...) base::message(sprintf(...))
  } else {
    function(...) {}
  }
}

##' Create a matrix from a list by rbinding all columns together
##' @title Create matrices from lists
##' @param x A list, or something coercable to a list
##' @export
rbind_list <- function(x) {
  do.call("rbind", as.list(x))
}
##' @export
##' @rdname rbind_list
cbind_list <- function(x) {
  do.call("cbind", as.list(x))
}

##' Spline interpolation in log-x space
##' @title Spline interpolation in log-x space
##' @param x,y Vectors giving coordinates of points to be
##' interpolated.  The x points should be naturally on a log scale,
##' and for \code{splinefun_loglog} both x and y should be on a log
##' scale.
##' @param ... Additional parameters passed to
##' \code{\link{splinefun}}.
##' @export
##' @author Rich FitzJohn
splinefun_log <- function(x, y, ...) {
  f <- splinefun(log(x), y, ...)
  function(x) {
    f(log(x))
  }
}
##' @export
##' @rdname splinefun_log
splinefun_loglog <- function(x, y, ...) {
  f <- splinefun(log(x), log(y), ...)
  function(x) {
    exp(f(log(x)))
  }
}

vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}
viapply <- function(X, FUN, ...) {
  vapply(X, FUN, integer(1), ...)
}
vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}
vcapply <- function(X, FUN, ...) {
  vapply(X, FUN, character(1), ...)
}