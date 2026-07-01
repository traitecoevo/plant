# Trait gradient of grow_individual_to_size (#472 scope B, the last FF16 surface).
# CI-runnable in plain R: ff16_grow_to_size_gradient_impl is compiled into plant.so
# with the XAD adjoint tape resolved at load against odelia, so no on-the-fly
# compilation is needed. The honest-scope gap to the fully-adaptive live solver (the
# frozen-grid caveat) is described in the AutoDiff guide.

test_that("grow_individual_to_size_gradient reconstructs grow_individual_to_size", {
  s    <- FF16_Strategy()
  indv <- Individual("FF16", "FF16_Env")(s)
  env  <- Environment("FF16")
  targets <- c(2, 5, 10)

  g   <- grow_individual_to_size_gradient(indv, targets, "height", env, time_max = 200)
  ref <- grow_individual_to_size(indv, targets, "height", env, time_max = 200)

  # The replay reproduces t* / state up to the live uniroot tolerance (R/individual.R
  # uniroot()s t* with default tol ~ 1.2e-4, so ~1e-4 in t / ~1e-4 in state).
  expect_equal(as.numeric(g$time), as.numeric(ref$time), tolerance = 5e-4)
  expect_equal(unname(g$state), unname(ref$state), tolerance = 5e-4)
  expect_equal(dim(g$d_state), c(length(targets), 5L, length(ff16_default_traits())))
})

test_that("d(height)/d(theta) is the IFT identity (~0: height is pinned to target)", {
  indv <- Individual("FF16", "FF16_Env")(FF16_Strategy())
  g <- grow_individual_to_size_gradient(indv, c(2, 5, 10), "height", Environment("FF16"),
                                        time_max = 200)
  # The target size cancels in d y_c/d theta = partial + rate*dt*/dtheta for c = size.
  expect_lt(max(abs(g$d_state[, "height", ])), 1e-6)
})

test_that("AD matches a frozen-schedule central FD (the exact contract)", {
  s    <- FF16_Strategy()
  indv <- Individual("FF16", "FF16_Env")(s)
  env  <- Environment("FF16")
  targets <- c(2, 8)
  traits  <- c("a_p1", "lma", "rho", "a_b1", "a_l1", "k_l")

  g  <- grow_individual_to_size_gradient(indv, targets, "height", env, traits = traits,
                                         time_max = 200)
  sh <- grow_individual_bracket(indv, targets, "height", env, time_max = 200)$time
  pp <- unlist(s$pars)

  # FD over the SAME frozen-schedule function the AD differentiates: perturb the trait,
  # recompute h0 (the trait-dependent seedling size), hold the schedule sh fixed. A
  # trait flows through the whole trajectory, so the optimal FD step spans orders of
  # magnitude across traits; use a per-quantity plateau picker over a step ladder (the
  # rung most self-consistent with its neighbour), as the guide prescribes.
  y0_of <- function(pp2) {
    s2 <- FF16_Strategy()
    for (nm in names(pp2)) s2$pars[[nm]] <- pp2[[nm]]
    i2 <- Individual("FF16", "FF16_Env")(s2)
    stats::setNames(i2$ode_state, i2$ode_names)
  }
  val <- function(pp2) {
    v <- ff16_grow_to_size_gradient_impl(pp2, env, y0_of(pp2), sh, as.numeric(targets),
                                         0L, traits, TRUE)
    list(time = v$time, state = v$state)
  }
  ladder  <- c(3e-3, 1e-3, 3e-4, 1e-4)
  plateau <- function(seq) seq[which.min(abs(diff(seq)))]
  nonh    <- setdiff(colnames(g$state), "height")

  for (tr in traits) {
    k   <- match(tr, traits)
    fds <- lapply(ladder, function(rh) {
      h  <- rh * abs(pp[[tr]])
      fl <- val(`[[<-`(pp, tr, pp[[tr]] - h))
      fu <- val(`[[<-`(pp, tr, pp[[tr]] + h))
      list(time = (fu$time - fl$time) / (2 * h), state = (fu$state - fl$state) / (2 * h))
    })
    for (i in seq_along(targets)) {
      fd_t <- plateau(vapply(fds, function(z) z$time[i], 0))
      expect_equal(unname(g$d_time[i, k]), unname(fd_t), tolerance = 2e-3,
                   info = sprintf("d(t*)/d(%s) @ size %g", tr, targets[i]))
      for (c in nonh) {
        fd_s <- plateau(vapply(fds, function(z) z$state[i, c], 0))
        # Mixed abs/rel: rel guards big derivatives, abs the near-zero components
        # (e.g. fecundity before maturity).
        expect_lt(abs(g$d_state[i, c, k] - fd_s), 1e-6 + 2e-3 * abs(fd_s))
      }
    }
  }
})

test_that("grow_individual_to_size_gradient rejects non-FF16 and bad inputs", {
  indv <- Individual("FF16", "FF16_Env")(FF16_Strategy())
  env  <- Environment("FF16")
  expect_error(grow_individual_to_size_gradient(indv, c(5, 2), "height", env),
               "sorted")
  expect_error(grow_individual_to_size_gradient(indv, 5, "not_a_state", env),
               "ODE state names")
})
