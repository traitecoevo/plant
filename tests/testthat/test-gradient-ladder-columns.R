# The column-by-column referee: the sweep against a forward tangent of the
# same trajectory, one column at a time.
#
# Split from the rest of the two-species checks because it is the expensive half -- a tangent column
# is 4.8 s and this asks for twenty-two of them -- and because the claim is
# different in kind. Those ask whether accumulation across cohorts and
# species is right; this asks whether each column agrees with an independent
# derivative of the same run, which is the contraction the whole rung rests on.

test_that("the replay reaches the run's own census", {
  # Before any derivative the reference reports means anything, it has to be a
  # derivative of the same function. The tangent run replays the recorded step
  # SIZES rather than the times -- a size differenced back out of two recorded
  # times is not the size that was taken, since fl(fl(t + h) - t) != h -- and
  # rather than a controller of its own, which would differentiate the
  # controller. Landing on the run's own census is what says the replay worked.
  #
  # It is also every trajectory check's floor: the reference reports derivatives
  # at the state it actually reached.
  shared <- ladder_shared("two_by_two")
  stand <- shared$stand
  ladder_require_regime(stand, "stand")
  result <- shared$gradient
  ladder_report_margin("the replay reaches the run's own census",
                       ladder_replay_floor(stand, result$value), 1e-9)

  # Non-vacuity of the seed: an unseeded replay must report no derivative, or a
  # zero tangent would be evidence of nothing about a column that reads zero.
  quiet <- ladder_trajectory_tangent(
    stand, numeric(length(colnames(result$gradient))))
  expect_true(all(quiet$tangent == 0))
})

test_that("the gradient contracts to a tangent of the same trajectory", {
  # The reference traverses both reductions, the stage recursion and the
  # introduction boundary, and none of the transposes under test are on its path.
  # Here the Jacobian can no longer be formed, so the contraction returns.
  shared <- ladder_shared("two_by_two")
  stand <- shared$stand
  result <- shared$gradient
  columns <- colnames(result$gradient)

  direction <- ladder_seeds(length(columns),
                            scale = ladder_block_scale(result$gradient[1, ]))
  reference <- ladder_trajectory_tangent(stand, direction)$tangent
  observed <- as.vector(result$gradient %*% direction)
  expect_length(reference, nrow(result$gradient))

  # Normalised by the largest single contribution rather than by their sum. The
  # sum cancels -- that cancellation is the contraction's blind spot, and it is
  # why the coordinate directions below are run beside it -- so dividing by it
  # would report cancellation as error.
  term <- max(abs(sweep(result$gradient, 2, direction, `*`)))
  ladder_report_margin("the contraction against the trajectory tangent",
                       max(abs(reference - observed)) / term,
                       3 * ladder_trajectory_agreement())
})

test_that("each shortlisted trait's own column is refereed", {
  # A contraction against one direction masks a column whose true value is near
  # zero, so every shortlisted parameter is asked for as a coordinate direction:
  # one replay, one exact column of the census-by-trait Jacobian.
  #
  # One bound for every column. There were two, the looser one held the pair of
  # allometric constants, and it was the class the pair belongs to that was wrong
  # rather than the pair -- so a second bound would now excuse the same shape again.
  shared <- ladder_shared("two_by_two")
  stand <- shared$stand
  result <- shared$gradient
  columns <- colnames(result$gradient)
  present <- columns[ladder_bare_traits(columns) %in% ladder_shortlist()]
  # Every birth-size parameter is a column here, which is the claim the count
  # asserts: a shortlist that silently lost one would read as a shorter loop.
  expect_true(all(ladder_birth_size_parameters() %in%
                    ladder_bare_traits(present)))
  expect_length(present, 22L)

  for (name in present) {
    reference <- ladder_trajectory_tangent(
      stand, ladder_trait_direction(columns, name))$tangent
    peak <- max(abs(result$gradient[, name]))
    expect_gt(peak, 0)
    ladder_report_margin(paste("column", name),
                         max(abs(reference - result$gradient[, name])) / peak,
                         ladder_trajectory_agreement())
  }
})
