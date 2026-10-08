# Pricing a declared zero, which needs a reference that does not declare it.
#
# Every other check in this ladder referees the sweep against a forward tangent.
# The tangent is exact and it settles whether the transpose is the transpose of its
# forward function -- but it inherits every equation the model imposes, the seed
# height's among them, so the two agree there for free and the agreement reads as a
# pass. A whole-run difference rebuilds the strategy, which re-runs preparation, so
# it carries the channel both differentiated paths hold still.
#
# This is the shortlist of columns that channel reaches, and no more: it costs two
# model runs per column.

test_that("a differenced column moves one registered parameter", {
  # The hyperparameter function derives registered parameters from traits, so a
  # difference taken on the trait vector answers what a species further along a
  # trait spectrum does while the sweep answers what one parameter does. Recorded
  # here as data rather than as prose, because which traits fan out is a property
  # of that function and changes when it does.
  traits <- ladder_traits()$fast
  for (name in names(traits)) {
    fanout <- ladder_trait_fanout(traits, name)
    message(sprintf("  %-5s also moves: %s", name,
                    if (length(fanout)) paste(names(fanout), collapse = ", ")
                    else "(nothing)"))
  }
  # lma is the one the fixture carries, and it is why the guard exists.
  expect_gt(length(ladder_trait_fanout(traits, "lma")), 0)
  expect_length(ladder_trait_fanout(traits, "a_l1"), 0)

  # The guard's sensitivity, established by breaking what it watches: no other
  # check in this suite notices a difference taken along a direction, because the
  # answer stays finite, plausible and correctly signed.
  before <- c(lma = 0.0825, k_l = 2.038, r_l = 572.4)
  expect_error(
    ladder_assert_one_parameter(before, c(lma = 0.0826, k_l = 2.030, r_l = 572.4),
                                "lma"),
    "moved 2 registered parameters")
  expect_silent(
    ladder_assert_one_parameter(before, c(lma = 0.0826, k_l = 2.038, r_l = 572.4),
                                "lma"))
})

test_that("the whole-run difference is in its own domain on this fixture", {
  # The check on the check. A re-run difference is unusable at production, so its
  # validity is established on this fixture rather than assumed.
  traits <- ladder_traits()$fast
  for (name in ladder_birth_size_parameters()) {
    got <- ladder_run_difference_stable(traits, name)
    message(sprintf(
      "\n  %-5s plateau at %s (step %.0e), spread %.2e over the pair and %.2e over all four; seed-height slope %+.4e",
      name, got$plateau_at, got$step, got$spread, got$spread_all,
      got$seed_height_slope))
    # ASSERTED, not skipped on. The fixture is written and the step ladder is
    # fixed, so the reference resolving here is a property of the pair -- and it
    # ceasing to is the model having moved under a fixture that did not.
    expect_lt(got$spread, 1e-3)
    # Non-vacuity: the channel this reference exists to price has to be live, or
    # agreeing with the sweep would prove nothing about the declaration.
    expect_gt(abs(got$seed_height_slope), 0)
  }
})

test_that("the birth-size channel is priced rather than asserted", {
  # Both differentiated paths impose d(seed height)/d(trait) = 0, and it is the one
  # term no other instrument in this suite can referee. This is that instrument.
  # What it reports is a ratio, not a pass: the declaration is a modelling choice,
  # and what a suite owes it is a number.
  #
  # The number moves with run length, so the run length is reported beside it --
  # a_l1 is a factor of two out at four tenths of a year and three parts in a
  # thousand at four years, and quoting either without the lifetime is quoting
  # nothing.
  traits <- ladder_traits()$fast
  stand <- ladder_stand_allometric_probe(TRUE)
  gradient <- stand_gradient(stand)
  cols <- colnames(gradient$gradient)

  message("\n  sweep against a whole-run difference, leaf area, 0.4 yr:")
  for (name in ladder_birth_size_parameters()) {
    got <- ladder_run_difference_stable(traits, name)
    if (got$spread > 1e-3) {
      message(sprintf(
        "    %-5s difference out of its own domain here (spread %.1e, plateau at %s)",
        name, got$spread, got$plateau_at))
      next
    }
    sweep <- gradient$gradient[1, paste0("1.", name)]
    ratio <- sweep / got$gradient[[1]]
    message(sprintf("    %-5s difference %+.6e  sweep %+.6e  ratio %6.3f",
                    name, got$gradient[[1]], sweep, ratio))
    expect_true(is.finite(ratio))
    # ⚠️ DELIBERATELY THREE ORDERS LOOSER THAN THE READING, and the looseness is
    # what makes it assertable. Every ratio reads 1.000 and 1e-04 would hold with
    # an order of margin, but the reference differences the same bisected
    # root-find the row check below bounds, and its floor is measured as a gap
    # between steps -- which reads the opposite of an error that rises as the
    # step falls. So a bound taken off that ladder would be a claim about a
    # convergence order this reference does not have. A bound of 1e-02 makes no
    # such claim and still rejects what could go wrong here: a dropped channel
    # (0 or infinite), a sign flip (-1), a row at a fixed fraction (0.5).
    expect_lt(abs(ratio - 1), 1e-2)
  }
})

test_that("the seed's geometry row is refereed against a rebuilt strategy", {
  # The birth height is not computed forwards -- it is the root of
  # mass_live_given_height(h) = omega, solved off the tape, and its derivative is
  # declared by the implicit function theorem. So the row has to be checked
  # against the condition it claims to solve, by rebuilding the strategy at a
  # perturbed parameter and finding the root again.
  #
  # The introductions file prices the same claim for two parameters through a whole rate
  # evaluation. This one has no patch, no environment and no rates on the path,
  # so a disagreement localises to the recording rather than to anything after it.
  patch <- ladder_patch_one()
  ladder_require_regime(patch, "patch")
  columns <- ladder_trait_names_tf24(patch)

  for (name in ladder_birth_size_parameters()) {
    at <- match(paste0("1.", name), columns)
    expect_false(is.na(at))
    # Named rather than defaulted: the bound below is the reference's own error at
    # THIS step, so the two have to be one number.
    rel <- 1e-4
    got <- ladder_seed_geometry_tangent_tf24(patch, at)
    ref <- ladder_boundary_difference(patch, name, rel = rel)

    # The recording returns the root itself, so the value is the rebuild's own.
    ladder_report_margin(
      sprintf("the seed height's value, %s", name),
      abs(got$height - ref$value[["height"]]) / abs(ref$value[["height"]]),
      1e-10)

    # Non-vacuity first: each of these reaches the residual, so a zero here is a
    # dropped channel rather than a small disagreement.
    expect_gt(abs(got$dheight), 0)

    # And the row, against a difference of the rebuilt root. What bounds it is
    # where the root-find stopped, not how fine the step is: the reference
    # differences a bisection whose answer is quantised at half a bracket, and a
    # relative step of 1e-4 in the parameter is a step of a few nanometres in the
    # height, which is smaller than the quantisation it is dividing.
    #
    # 1e-4 is also where this is the WHOLE bound: the reference's truncation there
    # is a ten-thousandth of it for every one of the eight. At 1e-3 truncation
    # returns at a seventh and the bound would need a second term.
    #
    # No factor of ten. This is a hard interval, not a scale estimate, so over it
    # means the row disagrees by more than a bisected root-find can account for.
    ladder_report_margin(
      sprintf("the seed height's row, %s", name),
      abs(got$dheight - ref$row[["height"]]) / abs(ref$row[["height"]]),
      ladder_bisection_bound(
        patch, 1L, ref$value[["height"]],
        abs(ladder_strategy_parameter(patch, 1L, name)) * rel,
        ref$row[["height"]]))
  }

  # And the structural half, over EVERY other column rather than a written pair.
  # The residual reads the birth-size parameters and no others, so the claim is
  # about the complement -- and a written pair is a claim about two names that
  # goes stale when the set moves. It went stale once already: the pair was
  # c("k_I", "a_l1") and `a_l1` is a birth-size parameter, so the loop asserted on
  # one name and said two.
  others <- setdiff(ladder_bare_traits(columns), ladder_birth_size_parameters())
  expect_gt(length(others), length(ladder_birth_size_parameters()))
  for (name in others) {
    at <- match(paste0("1.", name), columns)
    expect_identical(ladder_seed_geometry_tangent_tf24(patch, at)$dheight, 0,
                     info = paste("reaches the seed-height residual:", name))
  }
})

# The claim these lists make is about NUMBERS -- that certain columns come back
# exactly zero -- so that is what is checked. DO NOT check it against a label the
# C++ carries into every answer, which is a weaker subject: a label can agree with
# this list while the number disagrees with both, and it would have to be
# maintained beside the model that produces the number.
test_that("the columns declared zero are exactly zero, and no others are", {
  stand <- ladder_stand_two_by_two()
  g <- stand_gradient(stand)
  expect_false(any(stand_gradient_refused(g)))

  # Both directions. Every column the ladder declares zero is zero, and every
  # column that is zero is declared -- the second is what turns an unexplained
  # zero into a failure rather than a number nobody looks at.
  declared <- ladder_zero_by_construction()
  zero <- apply(g$gradient, 2, function(col) all(col == 0))
  observed <- ladder_bare_traits(names(zero)[zero])

  expect_setequal(observed, declared)
  message(sprintf("  %d of %d columns exactly zero: %s",
                  sum(zero), ncol(g$gradient),
                  paste(sort(unique(observed)), collapse = ", ")))

  # Non-vacuity: the zero columns really are a strict minority of a mostly-live
  # matrix, and the declared set is not the whole of it.
  expect_gt(sum(!zero), sum(zero))
  expect_true(all(is.finite(g$gradient)))
})
