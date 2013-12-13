list.to.array <- function(x) {
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

pad.list.to.array <- function(x)
  list.to.array(pad.matrix(x))

last <- function(x)
  x[[length(x)]]

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
run.cached <- function(expr, filename, regenerate=FALSE) {
  if (file.exists(filename) && !regenerate) {
    res <- readRDS(filename)
  } else {
    res <- eval.parent(substitute(expr))
    saveRDS(res, filename)
  }
  res
}
