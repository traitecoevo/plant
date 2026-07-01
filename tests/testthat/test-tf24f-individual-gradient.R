# TF24f individual grow-to-size trait gradient (#472 scope B, the "individuals" surface).
# A single plant grown in a FIXED environment to target sizes; the gradient is the
# finite-difference d(t*)/d(theta) and d(state at t*)/d(theta) over
# grow_individual_to_size (no resident feedback, no canopy/density -- the lightest
# surface). This is the prototype + the reference an exact-AD version must reproduce;
# the tracked collar (opt_root_psi_state) is re-evolved inside each grow, so the FD
# captures its theta-response automatically.

test_that("tf24f individual grow-to-size FD gradient is well-formed and sensible", {
  s <- TF24f_Strategy()
  ind <- TF24f_Individual(s)
  env <- Environment("TF24f")
  env$set_fixed_environment(1.0, 1e4)
  sizes <- c(2, 5, 9)
  g <- plant:::tf24f_grow_individual_to_size_gradient_fd(
    ind, sizes, "height", env, traits = c("vcmax_25", "lma", "K_s"), time_max = 400)

  # Shapes + the tracked collar is one of the returned state components.
  expect_length(g$time, length(sizes))
  expect_equal(dim(g$d_time), c(length(sizes), 3L))
  expect_true("opt_root_psi_state" %in% dimnames(g$d_state)[[2]])
  expect_true(all(is.finite(g$d_time)) && all(is.finite(g$d_state)))

  # Reconstructed stopping times match grow_individual_to_size.
  ref <- grow_individual_to_size(ind, sizes, "height", env, time_max = 400)
  expect_equal(g$time, ref$time, tolerance = 1e-6)

  # Physically-signed sensitivities of the time-to-size: costlier leaves (higher lma)
  # slow growth (d t* / d lma > 0); higher stem conductance (K_s) speeds it
  # (d t* / d K_s < 0). Checked at the largest target where the signal is clear.
  expect_gt(g$d_time[length(sizes), "lma"], 0)
  expect_lt(g$d_time[length(sizes), "K_s"], 0)
})

test_that("tf24f individual grow-to-size AD gradient matches the FD prototype (R1 tape)", {
  # The reverse-mode AD refine (build-order step 4): the 6-state grow {5 demog + tracked
  # collar} replayed over the frozen Cash-Karp schedule with the collar curvature-
  # linearised, plus the stopping-time IFT. It must reproduce the FD prototype -- but the
  # AD is the more accurate of the two: the collar relaxes from its clamped birth value
  # (0) over the early grow, a steep transient where the FD over the recon is noisy (its
  # d(t*)/d(lma) at the smallest target swings 1.35 -> 4.66 between 1e-5 and 1e-4 steps,
  # while AD gives a stable 4.67). So validate at the LARGEST target, where the collar has
  # equilibrated and the FD is converged; there AD matches FD to <1%.
  s <- TF24f_Strategy(); ind <- TF24f_Individual(s)
  env <- Environment("TF24f"); env$set_fixed_environment(1.0, 1e4)
  sizes <- c(2, 5, 9); tr <- c("vcmax_25", "lma", "K_s")
  ad <- plant:::tf24f_grow_individual_to_size_gradient_ad(ind, sizes, "height", env,
          traits = tr, time_max = 400)
  fd <- plant:::tf24f_grow_individual_to_size_gradient_fd(ind, sizes, "height", env,
          traits = tr, time_max = 400)

  # Reconstructed t* matches grow_individual_to_size exactly (the discovery pass is the
  # same bisection), and the AD output is well-formed.
  ref <- grow_individual_to_size(ind, sizes, "height", env, time_max = 400)
  expect_equal(as.numeric(ad$time), ref$time, tolerance = 1e-6)
  expect_equal(dim(ad$d_time), c(length(sizes), length(tr)))
  expect_true(all(is.finite(ad$d_time)) && all(is.finite(ad$d_state)))
  expect_true("opt_root_psi_state" %in% dimnames(ad$d_state)[[2]])

  # At the largest target the collar has equilibrated and the FD is converged: AD == FD.
  g <- length(sizes)
  expect_equal(ad$d_time[g, ], fd$d_time[g, ], tolerance = 0.01)
  # The well-behaved (non-cancelling) state-at-t* sensitivities match too.
  for (cc in c("mass_heartwood", "opt_root_psi_state"))
    expect_equal(ad$d_state[g, cc, ], fd$d_state[g, cc, ], tolerance = 0.01)

  # Physically-signed time-to-size: costlier leaves slow growth (d t*/d lma > 0), higher
  # stem conductance speeds it (d t*/d K_s < 0).
  expect_gt(ad$d_time[g, "lma"], 0)
  expect_lt(ad$d_time[g, "K_s"], 0)
})

test_that("grow_individual_to_size_gradient dispatches TF24f to the AD tape", {
  s <- TF24f_Strategy(); ind <- TF24f_Individual(s)
  env <- Environment("TF24f"); env$set_fixed_environment(1.0, 1e4)
  sizes <- c(2, 5, 9); tr <- c("vcmax_25", "lma", "K_s")
  disp <- grow_individual_to_size_gradient(ind, sizes, "height", env, traits = tr,
                                           time_max = 400)
  direct <- plant:::tf24f_grow_individual_to_size_gradient_ad(ind, sizes, "height", env,
              traits = tr, time_max = 400)
  expect_equal(disp$d_time, direct$d_time)
  expect_equal(disp$d_state, direct$d_state)
})

test_that("tf24f individual grow gradient is strategy-guarded", {
  ind <- FF16_Individual(FF16_Strategy())
  env <- FF16_Environment(); env$set_fixed_environment(1.0, 1e4)
  expect_error(
    plant:::tf24f_grow_individual_to_size_gradient_fd(ind, c(2, 5), "height", env),
    "TF24f strategy only")
})
