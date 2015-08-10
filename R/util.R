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
  d <- dim(x[[1]])
  if (is.null(d)) {
    if (length(unique(lapply(x, length))) > 1L) {
      stop("More than one length")
    }
    d <- length(x[[1]])
  }

  dimnames2 <- function(x) {
    dn <- dimnames(x)
    if (is.null(dn)) rep(list(NULL), length(dim(x))) else dn
  }
  array(unlist(x),
        c(d, length(x)),
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

##' Validate an object.  Currently only \code{Parameters} objects are
##' validated.
##' @title Validate an object
##' @param x Object
##' @param ... Additional arguments to be passed to methods
##' @export
validate <- function(x, ...) {
  UseMethod("validate")
}
##' @export
`validate.Parameters` <- function(x, ...) {
  plant <- parent.env(environment())
  ## TODO: This uses an implementation detail of RcppR6 that is not
  ## really OK to use; this could change at any moment.  Probably I'll
  ## expose this in some RcppR6 generated code eventually.
  type <- extract_RcppR6_template_type(x, "Parameters")
  get(sprintf("Parameters___%s__vdor", type), plant, inherits=FALSE)(x)
}

loop <- function(X, FUN, ..., parallel=FALSE) {
  if (parallel) {
    parallel::mclapply(X, FUN, ..., mc.preschedule=FALSE)
  } else {
    lapply(X, FUN, ...)
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

##' Make colours transparent
##' @title Make colours transparent
##' @param col Vector of colours
##' @param opacity Vector of opacities
##' @export
##' @examples
##' make_transparent("red", seq(0, 1, length.out=6))
##' make_transparent(c("red", "blue"), .5)
make_transparent <- function(col, opacity=.5) {
  alpha <- opacity
  if (length(alpha) > 1 && any(is.na(alpha))) {
    n <- max(length(col), length(alpha))
    alpha <- rep(alpha, length.out=n)
    col <- rep(col, length.out=n)
    ok <- !is.na(alpha)
    ret <- rep(NA, length(col))
    ret[ok] <- make_transparent(col[ok], alpha[ok])
    ret
  } else {
    tmp <- col2rgb(col)/255
    rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
  }
}

assert_named_if_not_empty <- function(x, name=deparse(substitute(x))) {
  if (length(x) > 0L) {
    nms <- names(x)
    if (is.null(nms) || any(nms == "") || any(duplicated(nms))) {
      stop(sprintf("All elements of %s must be uniquely named", name))
    }
  }
}

modify_list <- function(x, val) {
  modifyList(x, val[intersect(names(val), names(x))])
}

extract_RcppR6_template_type <- function(x, base) {
  cl <- class(x)[[1]]
  re <- sprintf("^%s<([^>]+)>$", base)
  if (!grepl(re, cl)) {
    stop("Unexpected type ", cl)
  }
  sub(re, "\\1", cl)
}
