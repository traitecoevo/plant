## Parallel version of Richardson extrapolation that will take
## advantage of tree's preference for computing multiple mutants at
## once.  Not yet set up to also do multiple x points at once.
##
## Defaults are set to match numDeriv::grad.
## Here's the simple version:
##
## gradient <- function(f, x,
##                      eps=1e-4, d=0.0001, r=4, v=2,
##                      zero_tol=sqrt(.Machine$double.eps/7e-7)) {
##   ## This is not actually hard to relax, I just don't need that right
##   ## now.  The big things are:
##   ##   - recycle d appopriately
##   ##   - x_up, x_down are matrices and are rbound/cbound together and
##   ##     we require that from f
##   ##   - once y is computed everything is the same.
##   if (length(x) != 1) {
##     stop("Requires scalar x")
##   }
##   ## Initial offset:
##   h0 <- abs(d * x) + eps * (abs(x) < zero_tol)
##   h <- h0 / (v^(seq_len(r) - 1))
##   x_up   <- x + h
##   x_down <- x - h
##
##   ## This is the bit I'm trying to get: *all* the x points are
##   ## computed at once.
##   y <- matrix(f(c(x_up, x_down)), r, 2)
##
##   ## Then bookkeeping galore.
##   a <- (y[,1] - y[,2]) / (2 * h)
##   m <- 1L
##   while (m < r) {
##     four_m <- 4 ^ m
##     a_next <- numeric(r - m)
##     for (i in seq_along(a_next)) {
##       a_next[i] <- (a[i + 1L] * four_m - a[i])/(four_m - 1.0)
##     }
##     a <- a_next
##     m <- m + 1L
##   }
##   a
## }

##' @export
gradient <- function(f, x,
                     eps=1e-4, d=0.0001, r=4, v=2,
                     zero_tol=sqrt(.Machine$double.eps/7e-7)) {
  ## This all corresponds to numDeriv's 'case 2'
  n <- length(x)

  ## Initial offset:
  h <- abs(d * x) + eps * (abs(x) < zero_tol)
  pts <- vector("list", r * n)
  dim(pts) <- c(r, n)
  dx <- matrix(NA, r, n)

  ## Note that this is not the most efficient way of computing these
  ## (note the do.call, rbind, and sapply constructs are all slow).
  ## But that does not matter in our case.
  j <- seq_len(n)
  for (k in seq_len(r)) {
    for (i in seq_len(n)) {
      dx_i <- h * (i == j)
      pts[[k, i]] <- rbind(x + dx_i, x - dx_i)
    }
    dx[k, ] <- 2 * h
    h <- h / v
  }

  xx <- do.call(rbind, pts)
  yy <- matrix(f(xx), 2L, r * n)
  a <- (yy[1,] - yy[2,]) / dx

  for (m in seq_len(r - 1L)) {
    four_m <- 4.0 ^ m
    a_next <- matrix(NA, r - m, n)
    for (i in seq_len(nrow(a_next))) {
      a_next[i,] <- (a[i + 1L,] * four_m - a[i,])/(four_m - 1.0)
    }
    a <- a_next
  }

  drop(a)
}

##' @export
bates_watts_D <- function(f, x,
                          eps=1e-4, d=0.0001, r=4,
                          zero_tol=sqrt(.Machine$double.eps/7e-7)) {
  v <- 2 # required
  len_f <- 1L # assume scalar output

  n <- length(x)  #  number of parameters (theta)
  h0 <- abs(d * x) + eps * (abs(x) < zero_tol)

  ## Points for computing first derivatives:
  pts <- vector("list", r * n)
  dx <- numeric(r * n)
  dim(pts) <- dim(dx) <- c(r, n)

  ## For second derivatives:
  pts2 <- dx2 <- vector("list", r * n * n)
  dim(pts2) <- dim(dx2) <- c(r, n, n)

  idx <- seq_len(n)
  h <- h0
  for (k in seq_len(r)) {
    for (i in seq_len(n)) {
      dx_i <- h * (i == idx)
      pts[[k, i]] <- rbind(x + dx_i, x - dx_i)
      for (j in seq_len(i - 1L)) {
        dx_j <- h * (j == idx)
        pts2[[k, j, i]] <- rbind(x + dx_i + dx_j,
                                 x - dx_i - dx_j)
        dx2[[k, j, i]] <- c(h[i], h[j])
      }
    }
    dx[k, ] <- h
    h <- h / v
  }

  xx <- rbind(x, do.call(rbind, pts), do.call(rbind, pts2),
              deparse.level=0)
  yy <- f(xx)

  ## Unpacking these points is really tricky, especially f2.
  i1 <- seq_len(2 * length(pts)) + 1L
  f0 <- yy[[1L]]
  f1 <- array(yy[i1], c(2L, dim(pts)))
  f2 <- yy[-c(1L, i1)]

  ## We have now done all the hard work to compute derivatives:
  D <- matrix(0, len_f, (n * (n + 3)) / 2)
  Daprox <- matrix(0, len_f, r)
  Hdiag  <- matrix(0, len_f, n)
  Haprox <- matrix(0, len_f, r)

  for (i in seq_len(n)) { # over parameters
    for (k in seq_len(r)) {
      Daprox[,k] <- (f1[1,k,i] -          f1[2,k,i]) / (2 * dx[k,i])
      Haprox[,k] <- (f1[1,k,i] - 2 * f0 + f1[2,k,i]) / dx[k,i]^2
    }
    for (m in seq_len(r - 1L)) {
      for (k in seq_len(r - m)) {
        Daprox[,k] <- (Daprox[, k + 1L] * (4^m) - Daprox[,k]) / (4^m - 1)
        Haprox[,k] <- (Haprox[, k + 1L] * (4^m) - Haprox[,k]) / (4^m - 1)
      }
    }
    D[,i]     <- Daprox[,1]
    Hdiag[,i] <- Haprox[,1]
  }

  u <- n
  idx <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      u <- u + 1L
      if (i == j) {
        D[,u] <- Hdiag[,i]
      } else {
        for (k in seq_len(r)) {
          tmp_f <- f2[c(idx, idx + 1L)]
          tmp_h <- dx2[[k, j, i]]
          Daprox[,k] <-
            (tmp_f[[1]] - 2 * f0 + tmp_f[[2]] -
             Hdiag[,i] * tmp_h[1L]^2 -
             Hdiag[,j] * tmp_h[2L]^2) /
               (2 * tmp_h[1L] * tmp_h[2L])
          idx <- idx + 2L
        }
        for (m in seq_len(r - 1L)) {
          for (k in seq_len(r - m)) {
            Daprox[,k] <- (Daprox[,k+1] * (4^m) - Daprox[,k]) / (4^m - 1)
          }
        }
        D[,u] <- Daprox[,1]
      }
    }
  }

  list(D=D, p=n, f0=f0, x=x)
}

## NOTE: I'm capitalising Hessian here to avoid colliding with numDeriv
##' @export
Hessian <- function(f, x, ...) {
  hessian_D(bates_watts_D(f, x, ...))
}

hessian_D <- function(obj) {
  D <- obj$D
  n <- obj$p
  H <- matrix(NA, n, n)
  H[upper.tri(H, TRUE)] <- D[-seq_len(n)]
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  H
}

##' @export
slope_info <- function(f, x, ...) {
  obj <- bates_watts_D(f, x, ...)
  list(x=x, fx=obj$f0, gr=obj$D[seq_along(x)], H=hessian_D(obj))
}
