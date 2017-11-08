##' Helper function for making bounds
##' @title Trait bounds
##' @param ... Named list, each element of which is a 2-element
##' numeric vector of lower and upper bounds.
##' @export
##' @examples
##' bounds(lma=c(0.01, 10))
##' bounds(lma=c(0.01, 10), rho=c(1, 1000))
bounds <- function(...) {
  x <- list(...)
  if (length(x) == 0) {
    stop("Need at least one argument")
  }
  if (!all(vapply(x, length, integer(1)) == 2)) {
    stop("All entries must be length 2")
  }
  if (is.null(names(x)) || any(names(x) == "")) {
    stop("All elements must be named")
  }
  ret <- rbind_list(x)
  colnames(ret) <- c("lower", "upper")
  ret
}

##' @param bounds A set of bounds
##' @param finite Logical indicating if bounds must be finite
##' @rdname bounds
##' @export
check_bounds <- function(bounds, finite=FALSE) {
  if (!is.matrix(bounds)) {
    stop("bounds must be a matrix")
  }
  if (ncol(bounds) != 2) {
    stop("bounds must have two columns")
  }
  if (is.null(rownames(bounds))) {
    stop("bounds must have rownames")
  }
  if (finite && any(!is.finite(bounds))) {
    stop("bounds must be finite")
  }
  colnames(bounds) <- c("lower", "upper")
  invisible(bounds)
}

##' @param trait_names Character vector of trait names
##' @rdname bounds
##' @export
bounds_infinite <- function(trait_names) {
  n <- length(trait_names)
  b <- cbind(lower=rep(-Inf, n), upper=rep(Inf, n))
  rownames(b) <- trait_names
  b
}

##' @export
##' @rdname bounds
##' @param x a point to detect if it lies within bounds
check_point <- function(x, bounds) {
  if (is.matrix(x)) {
    if (ncol(x) != nrow(bounds)) {
      stop("Invalid size x")
    }
  } else {
    if (length(x) != nrow(bounds)) {
      stop("Invalid size x")
    }
    x <- rbind(x, deparse.level=0)
  }
  if (is.null(names(x))) {
    colnames(x) <- rownames(bounds)
  } else if (names(x) != rownames(bounds)) {
    stop("Incorrect names on x")
  }
  tx <- t(x)
  if (any(tx < bounds[, "lower"] | tx > bounds[, "upper"])) {
    stop("Value does not lie within bounds")
  }
  invisible(x)
}

##' @importFrom stats uniroot
positive_1d <- function(f, x, dx, lower=-Inf, upper=Inf, tol=1e-3) {
  root <- function(b, type) {
    x <- b[[type]]$x
    fx <- b[[type]]$fx
    if (prod(fx[1:2]) < 0) {
      ## The suppressWarnings is here to stop the warning
      ##   -Inf replaced by maximally negative value
      ## which we're actually OK with.
      suppressWarnings(uniroot(f, x,
                               f.lower=fx[[1]], f.upper=fx[[2]],
                               tol=tol)$root)
    } else {
      if (type == "lower") x[[1]] else x[[2]]
    }
  }
  b <- positive_1d_bracket(f, x, dx, lower, upper)
  c(root(b, "lower"), root(b, "upper"))
}

positive_1d_bracket <- function(f, x, dx, lower, upper, grow=2) {
  fx <- f(x)
  if (fx < 0) {
    stop("Don't yet support doing this with no positive values")
  }

  L <- U <- x
  fL <- fU <- fx
  dL <- dU <- dx

  bracket <- function(x, dx, bound) {
    cleanup <- function(x, x_next, fx, fx_next) {
      if (dx < 0) {
        x <- c(x_next, x)
        fx <- c(fx_next, fx)
      } else {
        x <- c(x, x_next)
        fx <- c(fx, fx_next)
      }
      list(x=x, fx=fx)
    }
    hit_bounds <- FALSE
    repeat {
      x_next <- x + dx
      if ((dx < 0 && x_next < bound) || (dx > 0 && x_next > bound)) {
        x_next <- bound
        hit_bounds <- TRUE
      }
      fx_next <- f(x_next)
      if (fx_next < 0 || hit_bounds) {
        return(cleanup(x, x_next, fx, fx_next))
      } else {
        x <- x_next
        fx <- fx_next
        dx <- dx * grow
      }
    }
  }

  list(lower=bracket(x, -dx, lower),
       upper=bracket(x, dx, upper))
}

## This is a multidimensional version of positive_1d.  It's a hack for
## now.
positive_2d <- function(f, x, lower, upper, n_total=200) {
  lower <- rep1(lower, length(x))
  upper <- rep1(upper, length(x))

  requireNamespace("plant.ml")
  is_nonnegative <- function(y) y >= 0.0
  pts <- lapply(seq_along(lower), function(i) c(lower[[i]], upper[[i]]))
  m0 <- rbind(x, unname(as.matrix(do.call("expand.grid", pts))),
              deparse.level=0)
  plant.ml::delaunay_run_map(m0, f, is_nonnegative,
                             n_total=n_total, exploit=50)
}
