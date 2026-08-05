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
# ~3.7 MPa) -- brittleness worth its own issue, but not a place to anchor a
# regression test.

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
tf24_dry_start <- function(model = "TF24", lifetime = 10, theta = 0.005) {
  env <- Environment(model)
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(theta, 5))
  env$extrinsic_drivers_set_constant("rainfall", 1)

  p <- scm_base_parameters(model)
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  list(p = p, env = env)
}

## The completed dry-start run, computed once and shared. Two tests below need
## the same finished SCM -- the height ordering and the completion assertion --
## and each used to commission its own at ~2.9 s. Memoising is safe because both
## only read it; the pattern is the one used by short_run() in
## test-density-coordinate.R. Only the default (theta = 0.005) start is shared:
## the sweep below deliberately wants a fresh run per soil moisture, and must be
## free to throw.
tf24_dry_scm <- local({
  cache <- list()
  function(model = "TF24", lifetime = 10) {
    key <- paste(model, lifetime, sep = "/")
    if (is.null(cache[[key]])) {
      x <- tf24_dry_start(model, lifetime)
      scm <- SCM(model, environment_type(model))(x$p, x$env, Control())
      scm$collect <- TRUE
      scm$run()
      cache[[key]] <<- scm
    }
    cache[[key]]
  }
})

test_that("height_max() is the tallest cohort, not the first node (#571)", {
  # Under water limitation the top of the size distribution converges into a band
  # narrower than the refiner's finest spacing, and cohorts in it cross: TF24's
  # reserve-gated growth (#517) makes dh/dt depend on a cohort's own storage, so
  # two cohorts born moments apart into a rapidly wetting soil need not stay in
  # size order. That breaks the decreasing-height ordering that height_max() used
  # to exploit by returning nodes.front(), which then reported a height *below*
  # the tallest and only living cohort and truncated the light spline's domain.
  sp <- tf24_dry_scm()$patch$species[[1]]
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
  expect_true(is.finite(tf24_dry_scm()$offspring_production))

  # #571 recorded a non-monotone failure set across nine soil moistures, which
  # is what ruled out any single threshold as the cause. With the cause fixed
  # (#574) the wide sweep is regression ballast: six extra runs (~17 s) to
  # re-assert what the three values bracketing the residual-moisture floor
  # already cover. The full set -- 0.005, 0.008, 0.0099, 0.010, 0.0101, 0.012,
  # 0.015, 0.02, 0.03 -- is in #571 if this ever needs widening again. theta =
  # 0.005 is the shared run asserted above.
  for (theta in c(0.0099, 0.010, 0.0101)) {
    y <- tf24_dry_start(theta = theta)
    expect_no_error(run_scm(y$p, y$env, collect = FALSE))
  }
})

## REMOVED: "a failed light spline reports the patch state that caused it".
##
## It asserted that a resolution failure names the patch state behind it -- the
## "Patch state at that height" / "cohorts within" / "node heights are NOT
## decreasing" / "height_max() reports" lines, and the split between
## zero-density nodes scrambling the quadrature grid and cohorts genuinely
## crossing.
##
## Since #574 fixed the cause, the reproducer completes, so its
## `skip_if(is.na(msg), ...)` fired every run: it paid for a full 2.8 s SCM run
## and then asserted nothing. A test that can only ever skip is not a guard, it
## is a bill.
##
## The message-building code it covered is still live, and the assertions are in
## git history on this file if the failure ever returns. To resurrect it as a
## real test rather than a conditional one, drive the message builder from a
## synthetic patch state (crossed heights, some zero-density nodes) instead of
## trying to provoke a genuine spline failure -- that would assert
## unconditionally and cost nothing.
