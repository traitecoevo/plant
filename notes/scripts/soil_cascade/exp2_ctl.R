# odelia's exact controller (ode_control.hpp) with plant's defaults (src/control.cpp)
CTL <- list(tol_rel = 1e-4, tol_abs = 1e-4, a_y = 1.0, a_dydt = 0.0,
            h_min = 1e-6, h_max = 5, h_init = 1e-6)

errlevel <- function(y, dydt, h, ctl = CTL) {
  ctl$tol_rel * (ctl$a_y * abs(y) + ctl$a_dydt * abs(h * dydt)) + ctl$tol_abs
}

# NB: `std::max(a, b)` is `(a < b) ? b : a`, which does NOT propagate NaN the way
# R's max() does. With a = NaN it returns NaN, but on the next element a finite a
# returns a and *wipes* the NaN. Replicating that faithfully matters: it is why
# only a NaN surviving to the end of the loop ever reached the reject branch
# (odelia#52 / odelia PR#54). Using R's max() here would overstate the bug.
cxx_max <- function(a, b) if (isTRUE(a < b)) b else a

adjust <- function(h, ord, y, yerr, dydt, ctl = CTL) {
  rmax <- .Machine$double.xmin
  for (i in seq_along(y)) {
    D0 <- errlevel(y[i], dydt[i], h, ctl)
    rmax <- cxx_max(abs(yerr[i]) / abs(D0), rmax)
  }
  S <- 0.9
  if (isTRUE(rmax > 1.1)) {
    rr <- S / rmax^(1 / ord); if (rr < 0.2) rr <- 0.2
    hn <- h * rr; if (hn < ctl$h_min) hn <- ctl$h_min
    list(h = hn, shrank = hn < h)
  } else if (isTRUE(rmax < 0.5)) {
    rr <- S / rmax^(1 / (ord + 1)); if (rr > 5) rr <- 5
    hn <- h * rr; if (hn > ctl$h_max) hn <- ctl$h_max
    list(h = hn, shrank = FALSE)
  } else {
    list(h = h, shrank = FALSE)
  }
}

# Cash-Karp stage evaluation where uptake is recomputed from state at every
# stage (as plant does via the leaf model), so a bad probe can throw mid-step.
stage_with_uptake <- function(y, h, variant, uptake_fn, p = P) {
  b  <- list(c(), 1/5, c(3/40, 9/40), c(3/10, -9/10, 6/5),
             c(-11/54, 5/2, -70/27, 35/27),
             c(1631/55296, 175/512, 575/13824, 44275/110592, 253/4096))
  c5 <- c(37/378, 0, 250/621, 125/594, 0, 512/1771)
  c4 <- c(2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4)
  ff <- function(th) f(th, rain = 0, uptake = uptake_fn(th), variant, p)
  k <- vector("list", 6); k[[1]] <- ff(y)
  for (s in 2:6) {
    ys <- y; for (j in seq_len(s - 1)) ys <- ys + h * b[[s]][j] * k[[j]]
    k[[s]] <- ff(ys)
  }
  y5 <- y; y4 <- y
  for (s in 1:6) { y5 <- y5 + h * c5[s] * k[[s]]; y4 <- y4 + h * c4[s] * k[[s]] }
  list(y = y5, yerr = y5 - y4, dydt_out = k[[6]])
}
