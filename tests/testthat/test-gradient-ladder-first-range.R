# The range of the recording below the first introduction.
#
# Every other stand introduces at the initial time, so its first introduction sits at
# the first recorded step and the lowest range holds no step. A walk that ran
# that range and a walk that began above it produce the same numbers on every
# such fixture, which is why nothing else in this ladder constrains the lower end
# of the walk. This fixture resumes from a patch that already carries cohorts, so
# half its recording is below the first introduction, and three things here say the
# walk reached it: it swept one more range, it asked the inflow boundary for the
# steps in that range, and what it ended holding is the census's sensitivity to
# the first recorded state rather than to the state at the introduction.

test_that("the fixture holds recorded steps below its first introduction", {
  # Non-vacuity before anything else, and it is the property no other fixture in
  # this suite has.
  stand <- ladder_stand_resumed()
  ladder_require_regime(stand, "stand")
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  introduction <- which(diff(widths) != 0)
  message(sprintf("\n  %d recorded steps, %d introductions, widths %s",
                  length(widths), length(introduction),
                  paste(unique(widths), collapse = " -> ")))
  expect_gte(length(introduction), 2L)
  # Counted from one here, so an introduction at the first recorded step reads as 1.
  expect_gt(introduction[[1]], 1L)
  # And the first range is a substantial part of the run rather than a step or
  # two, so a walk that skipped it misses a term rather than a rounding.
  expect_gt(introduction[[1]], 0.25 * length(widths))

  # The state below the introduction carries cohorts, not the soil alone. On a run
  # from bare ground it would carry the environment only, and the census's
  # sensitivity to it is then twelve orders below its sensitivity one range up.
  environment_width <- stand$patch$ode_size - stand$patch$node_ode_size
  expect_gt(length(trajectory[[1]]$state), environment_width)
})

test_that("the sweep runs the range below the first introduction", {
  stand <- ladder_stand_resumed()
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  n_widening <- sum(diff(widths) != 0)
  steps <- ladder_step_count(trajectory)

  # One sweep, and both counts are its own return. DO NOT take one from a call and
  # the other off the solver afterwards: that pair can describe two different
  # sweeps.
  stand_gradient(stand)
  counts <- ladder_boundary_evaluations_tf24(stand)
  ranges <- counts$ranges

  # One range per width, and every width here has a step in it, so the count is
  # one MORE than the number of introductions. A walk beginning at the first introduction
  # reports exactly the number of introductions.
  message(sprintf("  %.0f ranges over %d introductions, %.0f boundary evaluations",
                  ranges, n_widening, counts$evaluations))
  expect_equal(ranges, n_widening + 1)

  # And the boundary is evaluated once per stage of every step swept. A walk that
  # skipped the lowest range would divide this by the steps in it, which on a
  # fixture that widens at its first step is nothing at all. It does not scale
  # with the metrics: a step is recorded once and the recording swept per metric.
  stages <- 6L
  expect_equal(counts$evaluations, stages * steps)
})

test_that("the adjoint the walk ends holding is the census's sensitivity to the first recorded state", {
  # The walk's lower end, measured against a forward tangent of the same
  # trajectory from the same state. The tangent is a forward replay at the
  # recorded step sizes and shares no code with the sweep, so a disagreement is
  # the walk's.
  stand <- ladder_stand_resumed()
  metrics <- census_metric_names_tf24()

  # Through stand_gradient so a sweep that refuses this stand skips by the
  # ladder's one gate and a sweep that BREAKS is re-raised. The adjoint itself
  # comes back with the gradient that ended holding it.
  stand_gradient(stand)
  lambda <- do.call(rbind, census_trait_gradient_tf24(stand)$at_first_state)
  rownames(lambda) <- metrics

  base <- ladder_range_base_state(stand, 0L)
  expect_equal(ncol(lambda), length(base))

  reference <- do.call(cbind, lapply(seq_along(base), function(i) {
    ladder_census_initial_state_tangent_tf24(
      stand, replace(numeric(length(base)), i, 1), 0L)$tangent
  }))
  rownames(reference) <- metrics

  # Non-vacuity: the state below the first introduction has to reach the census, or
  # the two matrices agree on a pair of zeros.
  expect_gt(max(abs(reference)), 1e-3)
  expect_gt(sum(reference != 0), length(base))

  ladder_report_margin("the walk's final adjoint against the forward tangent",
                       ladder_matrix_residual(lambda, reference),
                       ladder_trajectory_agreement())
})
