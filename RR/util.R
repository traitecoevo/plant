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
##' @keywords internal
validate <- function(x, ...) {
  UseMethod("validate")
}

##' @export
`validate.Parameters` <- function(x, ...) {
  plant <- parent.env(environment())
  ## TODO: This uses an implementation detail of RcppR6 that is not
  ## really OK to use; this could change at any moment.  Probably I'll
  ## expose this in some RcppR6 generated code eventually.
  types <- extract_RcppR6_template_types(x, "Parameters")
  constructor = do.call('sprintf', c("Parameters___%s__%s__vdor", types))
  get(constructor, plant, inherits=FALSE)(x)
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
##' @importFrom grDevices col2rgb rgb
##' @examples
##' util_colour_set_opacity("red", seq(0, 1, length.out=6))
##' util_colour_set_opacity(c("red", "blue"), .5)
util_colour_set_opacity <- function(col, opacity=.5) {
  alpha <- opacity
  if (length(alpha) > 1 && any(is.na(alpha))) {
    n <- max(length(col), length(alpha))
    alpha <- rep(alpha, length.out=n)
    col <- rep(col, length.out=n)
    ok <- !is.na(alpha)
    ret <- rep(NA, length(col))
    ret[ok] <- util_colour_set_opacity(col[ok], alpha[ok])
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

#' @importFrom utils modifyList
modify_list <- function(x, val) {
  modifyList(x, val[intersect(names(val), names(x))])
}

extract_RcppR6_template_types <- function(x, base) {
  cl <- class(x)[[1]]
  re <- sprintf("^%s<([^>]+)>$", base)
  if (!grepl(re, cl)) {
    stop("Unexpected type ", cl)
  }
  # Return a vector of type name strings
  as.list(strsplit(sub(re, "\\1", cl), ',')[[1]])
}

rep1 <- function(x, length.out, name=deparse(substitute(x))) {
  if (length(x) == length.out) {
    x
  } else if (length(x) == 1L) {
    rep_len(x, length.out)
  } else {
    stop(sprintf("%s must be length %d or scalar", name))
  }
}
