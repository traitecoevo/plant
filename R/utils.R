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
