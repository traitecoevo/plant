# Failures in the arid corner (#549 / #550 family), and honest reporting of them.
#
# Two separate defects, found while working out what plant needs to run a
# 416 mm dryland site (Dahra, Verbruggen et al. 2025):
#
#   1. TF24f's exact-gradient path entered the ci root-find in the hydraulic
#      shut-down state, where the solver throws rather than returning
#      non-finite, so a dry patch killed the whole run.
#   2. The adaptive interpolator reported "as refined as currently possible"
#      for a non-finite target, blaming resolution for a NaN.
#
# Both are diagnostic/robustness fixes: they change behaviour only in states
# that previously threw.
#
# The guard in (1) is tested through the SCM rather than on Leaf directly.
# Reaching the reversed-gradient state on a bare Leaf means driving it to a
# psi_soil at which the setup calls themselves fail first
# (find_root_psi throws "invalid f_ri" at psi_soil = 4 MPa, and
# find_psi_stem_from_psi_root leaves the transport spline's domain above
# ~3.7 MPa) -- brittleness worth its own issue, which became #576.
#
# #576 turned out NOT to be a too-narrow spline domain. TF24f's exact-gradient
# path asked for an acclimation gradient in the hydraulic shut-down state, where
# root_collar_psi_ is set to -root_psi_crit -- a collar potential at which soil
# uptake is *negative*, so the stem potential carrying it would have to be wetter
# than saturation and the transport inverse genuinely has no solution. The
# out-of-domain error was correct; the question was wrong. See the tests at the
# bottom of this file.

test_that("a dry TF24f patch no longer aborts on the ci root-find", {
  # The reproducer that found it: 5 layers held below the residual floor with
  # 1 m/yr rainfall. This asserts the *specific* failure is gone. The run can
  # still fail further downstream -- when growth stalls, cohorts pile up at the
  # introduction height and the light spline hits a genuine discontinuity (see
  # the interpolator test below and Verbruggen/capabilities.md) -- so this
  # deliberately checks the error is not the psi_stem_to_ci one rather than
  # asserting the run completes.
  env <- Environment("TF24f")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.005, 5))
  env$extrinsic_drivers_set_constant("rainfall", 1)

  p <- scm_base_parameters("TF24f")
  p$max_patch_lifetime <- 10
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))

  msg <- tryCatch({
    run_scm(p, env, collect = TRUE)
    NA_character_
  }, error = function(e) conditionMessage(e))

  expect_false(isTRUE(grepl("psi_stem_to_ci failed", msg, fixed = TRUE)))
  expect_false(isTRUE(grepl("do not bracket the root", msg, fixed = TRUE)))
})

test_that("adaptive interpolation names a non-finite target", {
  # check_err() compares against NaN, and every NaN comparison is false, so a
  # single NaN made its interval permanently unacceptable: refinement ran to
  # max_depth and then reported a resolution limit. The message sent debugging
  # in the wrong direction, so a non-finite value now says so.
  #
  # The light spline is the adaptive interpolator's only production caller, and
  # it is fed a C++ lambda, so drive it through the exposed test hook.
  expect_error(
    test_adaptive_interpolator(function(x) if (x > 0.5) NaN else x, 0, 1),
    "non-finite")

  # A genuinely unresolvable but finite target still reports resolution, and now
  # says what was exhausted.
  expect_error(
    test_adaptive_interpolator(function(x) as.numeric(x > 0.5), 0, 1),
    "as refined as currently possible")

  # A smooth target is unaffected.
  expect_silent(test_adaptive_interpolator(function(x) sin(x), 0, 1))
})

test_that("a resolution limit says where refinement stalled", {
  # Naming the x is the whole diagnosis: refinement stalls on a feature of the
  # target, so the location points straight at whatever put the feature there.
  # Without it the message says only that some feature somewhere is too narrow
  # (#571 took an afternoon to localise by hand).
  msg <- tryCatch(
    test_adaptive_interpolator(function(x) as.numeric(x > 0.5), 0, 1),
    error = conditionMessage)

  expect_match(msg, "at x = ")
  x <- as.numeric(sub(".*at x = ([0-9.e+-]+):.*", "\\1", msg))
  expect_equal(x, 0.5, tolerance = 1e-3)
  # And what the target does across the interval it could not resolve.
  expect_match(msg, "jumps from 0 to 1")
  expect_match(msg, "interval\\(s\\) still unresolved")
})

## Helper: the #571 dry start. Five layers held below the residual floor with
## 1 m/yr rainfall, which is a normal dryland initial condition rather than an
## edge case (Dahra is 416 mm MAP).
tf24_dry_start <- function(model = "TF24", lifetime = 10) {
  env <- Environment(model)
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.005, 5))
  env$extrinsic_drivers_set_constant("rainfall", 1)

  p <- scm_base_parameters(model)
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  list(p = p, env = env)
}

test_that("height_max() is the tallest cohort, not the first node (#571)", {
  # Under water limitation the top of the size distribution converges into a band
  # narrower than the refiner's finest spacing, and cohorts in it cross: TF24's
  # reserve-gated growth (#517) makes dh/dt depend on a cohort's own storage, so
  # two cohorts born moments apart into a rapidly wetting soil need not stay in
  # size order. That breaks the decreasing-height ordering that height_max() used
  # to exploit by returning nodes.front(), which then reported a height *below*
  # the tallest and only living cohort and truncated the light spline's domain.
  x <- tf24_dry_start()
  scm <- SCM("TF24", "TF24_Env")(x$p, x$env, Control())
  scm$run()

  sp <- scm$patch$species[[1]]
  h <- sp$heights

  # The premise of the test: this state really does violate the ordering. If a
  # future change makes the size distribution well-behaved, this stops being the
  # case that needs guarding and the assertions below become vacuous.
  skip_if(all(diff(h) <= 0), "node heights no longer invert here; see #571")

  expect_equal(sp$height_max, max(h))
  expect_gt(max(h), h[1])       # i.e. the front node is *not* the tallest
})

test_that("a dry-start TF24 run completes (#571)", {
  # Was: died in the light spline, because compute_competition()'s early exit
  # dropped every node past the first one below the query height -- including the
  # one cohort with appreciable density -- putting a fictitious step in the
  # competition profile that the refiner could not resolve.
  x <- tf24_dry_start()
  out <- run_scm(x$p, x$env, collect = TRUE)
  expect_true(is.finite(out$offspring_production))

  # The whole non-monotone failure set recorded in #571, which ruled out any
  # single threshold as the cause.
  for (theta in c(0.005, 0.008, 0.0099, 0.010, 0.0101, 0.012, 0.015, 0.02, 0.03)) {
    y <- tf24_dry_start()
    y$env$set_soil_water_state(rep(theta, 5))
    expect_no_error(run_scm(y$p, y$env, collect = FALSE))
  }
})

test_that("a failed light spline reports the patch state that caused it", {
  # The refiner can only report a narrow feature in a function; the patch knows
  # that the function is a light profile over a set of cohort heights. Kept as a
  # guard on the message: the reproducer no longer fails, so this asserts the
  # diagnosis is present *if* a resolution failure ever comes back.
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.005, 5))
  env$extrinsic_drivers_set_constant("rainfall", 1)

  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 10
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))

  msg <- tryCatch({
    run_scm(p, env, collect = TRUE)
    NA_character_
  }, error = conditionMessage)

  skip_if(is.na(msg), "dry-start TF24 now completes; see #571")

  expect_match(msg, "Patch state at that height")
  expect_match(msg, "cohorts within")
  # The ordering violation is the actionable part, so it must be reported rather
  # than left to be rediscovered.
  expect_match(msg, "node heights are NOT decreasing")
  expect_match(msg, "height_max\\(\\) reports")
  # And it must separate the two very different readings of that violation:
  # zero-density nodes scrambling the quadrature grid (bookkeeping) versus live
  # cohorts crossing (which would mean the characteristics themselves crossed,
  # and the method of characteristics forbids that). Measured here it is the
  # former, but assert that both counts are reported rather than pinning which.
  expect_match(msg, "between two cohorts of non-zero density")
  expect_match(msg, "nodes have zero density")
})

## ---------------------------------------------------------------------------
## #576: the transport spline's domain, and who was asking
## ---------------------------------------------------------------------------

## The #576 reproducer: the dry start above, but at the narrow rainfall window
## where TF24f died and TF24 did not. Rainfall 0 and 0.06 both ran, so this is a
## knife-edge in the middle rather than "drier fails".
tf24_dry_window <- function(model = "TF24f", rain = 0.05, lifetime = 10) {
  env <- Environment(model)
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.005, 5))
  env$extrinsic_drivers_set_constant("rainfall", rain)

  p <- scm_base_parameters(model)
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  list(p = p, env = env)
}

test_that("TF24f completes in the dry rainfall window (#576)", {
  # Both values failed with "Extrapolation disabled and evaluation point outside
  # of interpolated domain", from psi_from_transpiration's LOWER end. The
  # neighbours are included to keep the window itself under test: if a future
  # change moves the knife-edge rather than removing it, the run that starts
  # failing should be one of these.
  for (rain in c(0.03, 0.04, 0.05, 0.06)) {
    x <- tf24_dry_window(rain = rain)
    expect_no_error(run_scm(x$p, x$env, collect = FALSE))
  }
})

test_that("both TF24f gradient methods survive hydraulic shutdown (#576)", {
  # The cause was an asymmetry, not a bad number. The finite-difference branch of
  # TF24f_Strategy::solve_leaf checked prepare_collar_solve's return value and
  # returned a zero gradient when the operating point had been forced by
  # feasibility handling; the AD branch called evaluate_root_collar_psi (which
  # hides that return value) and then asked for a gradient regardless. So the
  # test that matters is that the two branches agree here, not that one number
  # came out right.
  for (ad in c(TRUE, FALSE)) {
    x <- tf24_dry_window()
    s <- x$p$strategies[[1]]
    s$use_ad_gradient <- ad
    x$p$strategies[[1]] <- s
    expect_no_error(run_scm(x$p, x$env, collect = FALSE))
  }
})

test_that("an out-of-domain transport lookup names the spline (#576)", {
  # Naming the spline, the point and the domain is the whole diagnosis: odelia's
  # own message is a bare sentence with none of the three, and #576 was localised
  # by bisecting the four candidate call sites by hand. Driven on a bare Leaf
  # because the SCM reproducer no longer reaches the state.
  l <- Leaf(vcmax_25 = 96, jmax_25 = 157.44, c = 2.04, b = 3.457268,
            psi_crit = 5.91988, root_c = 2.680147, root_b = 3.898245,
            root_psi_crit = 5.870283, beta2 = 1.5, a = 0.3,
            curv_fact_elec_trans = 0.7, curv_fact_colim = 0.99,
            GSS_tol_abs = 1e-3, vulnerability_curve_ncontrol = 100,
            ci_abs_tol = 1e-3, ci_niter = 1000, g1_TF24 = 7.5,
            beta_R_H = 3.4e2, beta_R_V = 9.4e3)

  msg <- tryCatch(l$transpiration(50, 0), error = conditionMessage)
  expect_match(msg, "transpiration_from_psi", fixed = TRUE)
  expect_match(msg, "psi_stem = 50", fixed = TRUE)
  expect_match(msg, "beyond the upper end", fixed = TRUE)
  # The domain, so the reader can see whether widening it would even help.
  expect_match(msg, "\\[0, [0-9.]+\\]")
  # And which call was asking, since the same spline is read from four places.
  expect_match(msg, "in Leaf::transpiration", fixed = TRUE)

  # A non-finite point must NOT be turned into a domain complaint. The guard uses
  # the same comparison odelia does (`u < lo || u > hi`) rather than negating an
  # in-range test, so NaN falls through to the spline as it always has and comes
  # back as NA -- caught downstream by the isfinite(profit_) guards, which report
  # the whole collar solve rather than one lookup. Negating an in-range test here
  # would read as a tightening and is really a behaviour change: it breaks the
  # documented profit_psi_stem_TF(NA, .) -> NA contract in test-leaf.r.
  expect_true(is.na(l$transpiration(NaN, 0)))
})
