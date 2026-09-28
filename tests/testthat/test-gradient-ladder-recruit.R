# The inflow boundary measured rather than described: what a marginal recruit's
# gradient does in the limit, and how many times the boundary's own row enters.
#
# Split from the introductions checks because both drive their own fixtures. The recruit
# checks sweep one stand per establishment scale and the count check sweeps two run
# lengths, so neither can take the shared sweep the other rung-5 checks read, and
# keeping them here leaves that file free to share one.

test_that("a marginal recruit sits in the stiff band, which is what makes the rest count", {
  # The fixture assertion is derived rather than chosen. The establishment
  # probability's derivative peaks at a known production, and a recruit anywhere
  # else is in a flat region where the limit below would pass on the signal being
  # small rather than on the mathematics being right.
  #
  # At the shipped establishment constant this fixture's recruits sit at 36 and 21
  # times that production, with an establishment probability of 0.998 -- a
  # comfortable seedling, and the wrong place to ask the question.
  stand <- ladder_shared("marginal_recruit")$stand
  stiffness <- ladder_recruit_stiffness(stand)
  establishment <- ladder_recruit_establishment(stand)
  message(sprintf("\n  recruit at %.3f of the peak production, establishment %.5f",
                  stiffness, establishment))
  expect_true(all(stiffness > 1 / 3 & stiffness < 3))
})

test_that("a marginal recruit's gradient exists, is finite, and tends to zero", {
  # Both the cumulative mortality and the log density diverge logarithmically as
  # establishment goes to zero, and the sensitivity of the log density diverges
  # with them. But the log density reaches everything downstream only through the
  # stem number, whose own sensitivity tends to zero, and the census carries stem
  # number linearly. So the finite answer is reached by seeding the quantity that
  # has one, and a recruit that cannot pay for itself contributes at second order
  # -- exactly zero is the correct first-order answer.
  #
  # The probe is the establishment decay, which reaches the census through the
  # boundary condition and nothing else, so its column is the recruit's own
  # contribution rather than a stand-wide one.
  scales <- c(30, 60, 120, 400, 2000)
  establishment <- numeric(length(scales))
  decay <- numeric(length(scales))
  for (k in seq_along(scales)) {
    stand <- ladder_stand_marginal_recruit(scales[[k]])
    result <- stand_gradient(stand)
    expect_true(all(is.finite(result$gradient)))
    establishment[[k]] <- ladder_recruit_establishment(stand)[[1]]
    columns <- colnames(result$gradient)
    decay[[k]] <- max(abs(result$gradient[
      , ladder_bare_traits(columns) == "recruitment_decay"]))
    message(sprintf("  establishment %.6f   |recruitment_decay| %.4e",
                    establishment[[k]], decay[[k]]))
  }

  # Non-vacuity: the sweep has to actually approach the limit, or "tends to zero"
  # is a statement about one point.
  expect_lt(min(establishment), 1e-3)
  expect_gt(max(establishment) / min(establishment), 1e3)

  # It tends to zero, and it does so from above at every step rather than only
  # between the endpoints.
  expect_true(all(diff(decay) < 0))
  expect_lt(decay[[length(decay)]] / decay[[1]], 1e-2)
})

test_that("the boundary's own term enters once per stage of every step", {
  # Confirming a derivative is not confirming its use. Every other check on this
  # channel reads the boundary condition's row -- against a rebuild, against a
  # tangent, entry by entry through the introduction map -- and none of them reads
  # how many times that row is multiplied in. A row correct per evaluation and wrong
  # in what it is multiplied by is a different failure from a wrong row, and neither
  # differentiated path can see it because both apply the same multiplier.
  #
  # The distinction is between a channel that acts once per plant and a channel that
  # acts once per step. The boundary node stands at the seed's height for a whole
  # run, so its VALUE changes once per plant; its condition is re-evaluated at every
  # stage of every step, so its ROW enters that many times.
  # Only the second is a count, and this is the count.
  stages <- 6L   # the Cash-Karp pair's stage count, two of whose output weights vanish

  measure <- function(lifetime) {
    p <- ladder_parameters("fast", lifetime = lifetime)
    p$node_schedule_times <- list(0)
    stand <- ladder_run(p)
    trajectory <- stand$store_trajectory()
    widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
    got <- ladder_boundary_evaluations_tf24(stand)
    list(steps = ladder_step_count(trajectory),
         introductions = sum(diff(widths) != 0),
         metrics = got$metrics, evaluations = got$evaluations,
         placements = got$placements)
  }

  short <- measure(0.4)
  long <- measure(0.8)
  for (m in list(short, long)) {
    message(sprintf(
      "  %3d steps, %d introduction(s), %d metrics: %.0f evaluations, %.0f placements",
      m$steps, m$introductions, m$metrics, m$evaluations, m$placements))
    expect_gt(m$steps, 1L)
    expect_gt(m$metrics, 1L)

    # Once per stage of every step, and NOT once per metric: a step is recorded
    # once and the recording swept per metric, so the model is evaluated as often
    # for three metrics as for one. A boundary term moved to act once per
    # introduction would collapse this from hundreds to a handful; a term applied
    # once per step rather than once per stage would divide it by six; and a
    # recording taken per metric instead of per step would multiply it by the
    # metric count.
    expect_equal(m$evaluations, stages * m$steps)

    # And the record ENGAGED. The sweep re-runs the model over the states the run
    # already visited, so every operating point it needs was found once already and
    # is placed rather than searched for -- but a record that quietly places nothing
    # searches again and returns the same numbers, because the search answers from
    # the same state. This count is the only thing that separates the two, and
    # nothing was reading it.
    #
    # A floor rather than an equality: the sweep addresses five stages of every step
    # and each places one point per cohort, so this is expected to be several times
    # the step count and is asserted only to be past it.
    expect_gt(m$placements, m$steps)
  }

  # Non-vacuity: the count has to be a count. A constant would satisfy every
  # equality above on one fixture and say nothing.
  expect_gt(long$steps, short$steps)
  expect_gt(long$evaluations, short$evaluations)
  expect_gt(long$placements, short$placements)
  message(sprintf("  the count tracks the step count: %.0f over %d steps against %.0f over %d",
                  short$evaluations, short$steps, long$evaluations,
                  long$steps))
})
