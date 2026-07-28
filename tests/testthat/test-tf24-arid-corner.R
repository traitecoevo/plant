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

test_that("a failed light spline reports the patch state that caused it", {
  # The #571 reproducer. The refiner can only report a narrow feature in a
  # function; the patch knows that the function is a light profile over a set of
  # cohort heights, and that here the node list has lost its decreasing-height
  # ordering -- which is what puts the fictitious step in the profile, because
  # Species::compute_competition() early-exits on the first node below the query
  # height and so drops the taller, and only living, cohort beyond it.
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
