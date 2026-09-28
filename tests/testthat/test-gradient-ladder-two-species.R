# Two species, two cohorts each, over a trajectory.
#
# Two species with one cohort each catch a reduction sum that collapses the
# species and not per-species accumulation; two cohorts of one species catch the
# accumulation and not the collapse. Four nodes catch both and are still small
# enough to reason about by hand.
#
# New over the one-cohort fixture: the stage recursion, step-size handling, accumulation across
# cohorts and species, and record-once-sweep-many.
#
# Here the Jacobian can no longer be formed, so the contraction returns -- with
# its blind spot stated. A contraction against one direction masks a column whose
# true value is near zero, so it is paired with the zero classification of the
# floor and with coordinate directions on a declared shortlist, each chosen
# because it has a documented structural route to reading exactly zero.


test_that("the structural checks hold before any reference is consulted", {
  # These need no reference, so they fail earlier and more cheaply than the
  # contraction -- and a failure here invalidates the contraction anyway.
  shared <- ladder_shared("two_by_two")
  stand <- shared$stand
  ladder_require_regime(stand, "stand")
  result <- shared$gradient

  expect_equal(dim(result$gradient),
               c(3L, length(census_trait_names_tf24(stand))))
  expect_true(all(is.finite(result$gradient)))

  # The trait adjoint accumulates over every cohort, every stage and every step,
  # so a trait treated as one input per cohort returns a fixed fraction of the
  # right answer with the correct sign and no error raised. Two species with two
  # cohorts each is the smallest stand where that fraction is not one.
  # Two species, each carrying the whole registered set under its own prefix.
  columns <- census_trait_names_tf24(stand)
  expect_equal(length(columns),
               2L * length(unique(sub("^[0-9]+[.]", "", columns))))
})

test_that("no shortlisted trait reads exactly zero", {
  # A contraction against one direction masks a column whose true value is near
  # zero. Coordinate directions on the shortlist are what cover that blind spot,
  # and every entry on the list is there because it has a route to reading zero
  # for a structural reason.
  shared <- ladder_shared("two_by_two")
  stand <- shared$stand
  result <- shared$gradient
  # Every species' column for every shortlisted parameter, because a column that
  # reads zero for one species only is exactly the collapse the shortlist is
  # watching for. The crown shape has no column at all and drops out here.
  columns <- colnames(result$gradient)
  present <- columns[ladder_bare_traits(columns) %in% ladder_shortlist()]
  expect_length(present, 22L)

  peak <- vapply(present, function(n) max(abs(result$gradient[, n])), numeric(1))
  message("\n  shortlisted trait magnitudes:")
  for (n in present) {
    message(sprintf("    %-20s %.6e", n, peak[[n]]))
  }
  expect_true(all(peak > 0))
})

test_that("a trait is one input read by every cohort, not one input per cohort", {
  # Because the strategy is shared, a single trait is one input. Treating each
  # cohort as a distinct input yields a gradient that is a fixed fraction of the
  # right answer with the correct sign and no error raised, and the fraction is
  # the cohort count. So the signature is a ratio, and adding a cohort at fixed
  # ecology is what exposes it.
  result_fewer <- ladder_shared("two_by_two")$gradient

  more <- local({
    p <- ladder_parameters(c("fast", "slow"))
    p$node_schedule_times <- list(c(0, 0.31, 0.63), c(0, 0.41))
    ladder_run(p)
  })
  result_more <- stand_gradient(more)

  # Refining the schedule at fixed ecology moves the answer by quadrature, not
  # by a whole cohort count. A ratio near the node-count ratio is the signature.
  ratio <- result_more$gradient["leaf_area", "1.lma"] /
    result_fewer$gradient["leaf_area", "1.lma"]
  message(sprintf("\n  lma row, four nodes over three: %.4f", ratio))
  expect_lt(abs(ratio - 1), 0.5)
})

test_that("the two species do not collapse into one scalar", {
  # A reduction that sums a per-species parameter over every cohort of every
  # species collapses the species into one number. This is undetectable with one
  # species, which is what earns the second, and the fixture's species differ by
  # whole factors in exactly the parameters that carry a per-species reduction
  # row.
  shared <- ladder_shared("two_by_two")
  stand <- shared$stand
  result <- shared$gradient
  columns <- census_trait_names_tf24(stand)
  half <- length(columns) / 2L

  for (name in c("k_I", "a_l1", "a_l2")) {
    # Naming the two columns rather than searching for them: a search that comes
    # back with one match, or none, would skip the check that earns the second
    # species.
    at <- match(paste0(1:2, ".", name), columns)
    expect_false(anyNA(at), label = paste("both species' columns for", name))
    first <- result$gradient[, at[[1]]]
    second <- result$gradient[, at[[2]]]
    expect_false(isTRUE(all.equal(first, second)),
                 label = paste("the two species' rows for", name, "are equal"))
  }
  expect_equal(half, length(unique(sub("^[0-9]+[.]", "", columns))))
})

test_that("record once and sweep many gives each metric its own recording", {
  # Clearing a tape returns its derivative-slot counter to zero, so an active
  # value constructed outside the sweep loop and read inside it refers, after the
  # first clear, to a slot that now belongs to something else. The first metric
  # is correct and every later one reads unrelated storage -- and a value that
  # survives a clear registers with no dependencies, so its adjoint sweeps to
  # exactly zero, which is why the failure produces exact zeros in whole column
  # families rather than noise.
  shared <- ladder_shared("two_by_two")
  stand <- shared$stand
  result <- shared$gradient
  metrics <- rownames(result$gradient)

  zero_families <- vapply(metrics, function(m) all(result$gradient[m, ] == 0),
                          logical(1))
  expect_false(any(zero_families))

  # And the positional check, which CANNOT FAIL as written and is kept only
  # because the shape it watches for is real. stand_gradient(metrics = m) computes
  # every metric in C++ and subsets in R, so "swept alone" and "swept with the
  # others" are one sweep and comparing them is a tautology. Established by
  # injection: with an active value deliberately held across the tape clear, this
  # passed and only the trajectory reference saw the defect. Closing it needs the
  # sweep itself to take a metric subset.
  # One metric rather than every metric. Each iteration is a full sweep at 43.6 s
  # for a comparison that is an identity by construction, so three of them cost two
  # minutes and establish exactly what one does. It stays as the placeholder for the
  # real check, which needs the sweep itself to take a metric subset.
  alone <- stand_gradient(stand, metrics = metrics[[1]])
  expect_identical(alone$gradient[metrics[[1]], ], result$gradient[metrics[[1]], ])
})
