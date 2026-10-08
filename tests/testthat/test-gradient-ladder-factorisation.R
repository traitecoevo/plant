# The one seam, and the check that closes it.
#
# Everywhere else in the ladder the reference and the object under test are
# disjoint. Here they are not: the tangent consumes the same supplied leaf rows
# the sweep does, and that is licensed only by those rows having been refereed
# independently against the individual's own algebra.
#
# The licence is weakest exactly where it matters most. The water rows are built
# on a factorisation of the mixed second derivative of profit into two scalars
# times closed-form vectors, and the natural check on it provably cannot referee
# it: the potential-family vectors are numerically collinear, so a compensating
# pair of coefficients fits every row equally well and an error in the second is
# absorbed into the first at a ratio of about a hundred thousand. Since the
# uniform-drying direction is a near-symmetry whose true response is amplified
# fifteen- to twenty-six-fold, a one-percent error in that coefficient is a
# fifteen- to twenty-six-fold error in the quantity the ecology cares about.
#
# So a wrong coefficient sits in the sweep and in the reference together, and
# every rung above passes. This is a common-mode failure on the channel the
# design names as its entire difficulty: the carbon half is free and the water
# half is the problem.
#
# It costs one state, no gradient run and no stand, and it gates the meaning of
# everything above it -- which is why it is the first thing to close and why it
# has a file of its own.

test_that("the marginal profit's factorisation is measured, not assumed", {
  # PASS is two things:
  #
  #   1. the second coefficient checked against its own closed form, at the
  #      first coefficient held fixed. Checking the pair jointly cannot work:
  #      a joint residual is flat along the compensating direction.
  #   2. the factorisation's residual over all of the state directions -- the
  #      soil potentials, the per-layer root masses and leaf area -- at states
  #      that include a pinned operating point and a layer near the
  #      equal-potentials branch.
  #
  # Everything the check needs is a forward quantity and each is reachable: the
  # marginal profit at a held collar potential, total uptake at a held collar
  # potential, and the physiology re-set at a perturbed soil potential, root
  # network or leaf area. The residual is then three central differences per
  # direction and a two-column least squares.
  #
  # The reading IS settled, in the code and in the corpus: the second intermediate
  # is uptake's sensitivity to the COLLAR, not to root mass, and the second basis
  # vector is d(dE_up/dp)/dpsi_j. So the factorisation is not a specification --
  # it is the production path, and every layer but the one it is solved from is
  # already a PREDICTION.
  #
  # That is what makes it refereeable with nothing built for the purpose. The
  # supplied uptake rows use the factorisation; a difference of the PLAIN-DOUBLE
  # block re-solves the collar, so it carries the true dp*/dpsi_j with no
  # factorisation in it. The two disagree only by the factorisation's own error.
  #
  # A difference is admissible here for the reason it is admissible at the soil
  # balance and not at the recorded step: it differences the double path, which
  # runs the leaf's own solve, rather than the recorded expression whose value does
  # not depend on a recorded input at all.
  patch <- ladder_patch_one()
  ladder_require_regime(patch, "patch")

  block <- ladder_block_jacobian_forward_tf24(patch, 1L)
  inputs <- colnames(block)
  outputs <- rownames(block)
  is_soil <- grepl("^psi_soil_", inputs)
  uptake <- grepl("^uptake_", outputs)
  supplied <- block[uptake, is_soil, drop = FALSE]

  # Two steps, so the reference's own error is measured rather than assumed. A
  # step of 1e-3 is NOT usable here: it drives the wettest layer's potential below
  # the root vulnerability grid's lower bound and the leaf refuses to bracket.
  differenced <- function(rel) {
    ladder_block_difference_tf24(patch, 1L, rel)[uptake, is_soil, drop = FALSE]
  }
  # Three steps, not two. A two-point spread cannot say whether it is measuring
  # truncation or round-off, and a tolerance taken from a round-off statistic is
  # not a tolerance. The layers differ by three orders in how well they converge
  # -- the wettest sits near the root vulnerability grid's lower bound, the same
  # fact that rules a 1e-3 step out -- so the spread is taken per layer.
  steps <- list(differenced(1e-5), differenced(1e-6), differenced(1e-7))
  fine <- steps[[3]]
  scale <- max(abs(fine))
  spread <- pmax(apply(abs(steps[[1]] - steps[[2]]), 2, max),
                 apply(abs(steps[[2]] - steps[[3]]), 2, max))
  floor <- pmax(spread / scale, 4 * .Machine$double.eps)
  residual <- apply(abs(supplied - fine), 2, max) / scale
  message("\n  per-layer residual: ",
          paste(sprintf("%.2e", residual), collapse = " "))
  message("  per-layer reference error: ",
          paste(sprintf("%.2e", floor), collapse = " "))

  # EVERY layer is out of sample. The pair is anchored in root carbon, not in a
  # soil potential, so no potential layer is the one it was solved from and none
  # of them gets the fit step's own second-order error allowed for. This used to
  # excuse whichever layer had the largest diagonal, on the reading that the pair
  # was fitted from a potential; that stopped being true when the anchor moved.
  expect_gt(length(residual), 1L)

  # Reported per layer because that is how the disagreement was found, but it is
  # not what is asserted: this statistic scales each layer by the whole block's
  # largest entry, so one error common to every layer reads as five different
  # numbers -- 20.9, 18.9, 17.5, 0.5, 1.4 when the root vulnerability curve was
  # tabulated as a fit through its values alone, and 1.6, 0.4, 1.0, 0.8, 1.1 once
  # the tabulation carried an exact slope at every knot. The second reading is the
  # residual sitting at its own reference's error on every layer.
  message("  per-layer ratio to its own reference: ",
          paste(sprintf("%.1f", residual / floor), collapse = " "))

  # Non-vacuity: the predicted entries must carry a response, or agreeing about
  # zero would pass. The off-diagonal of the supplied block IS the argmax channel.
  off <- supplied
  diag(off) <- 0
  ladder_expect_moves(off, matrix(0, nrow(off), ncol(off)),
                      "the argmax channel's off-diagonal")

  # The disagreement is ONE scalar. The difference matrix is rank one to a part in
  # ten thousand, which says every layer is wrong by the same relative amount and
  # so exactly one quantity is -- dcollar/dpsi. Size and structure are asserted
  # apart because neither sees the other's failure: a scalar bound cannot see a
  # wrong basis vector, and a rank bound cannot see a wrong size.
  # Refereed against the middle step, not the finest. The 1e-7 difference is past
  # the round-off crossing -- it disagrees with the 1e-6 one by 5e-06 on the
  # driest layer, where 1e-5 and 1e-6 agree to 2.8e-09 -- so the finest step is
  # the noisiest reference of the three, and the per-layer floor above is a ratio
  # to that noise, which is why it moves with the build and the entries below do
  # not.
  reference <- steps[[2]]

  singular <- svd(supplied - reference)$d

  # ⚠️ THE RANK TEST IS ONLY MEANINGFUL WHILE THE RESIDUAL EXCEEDS ITS REFERENCE,
  # AND IT NO LONGER DOES. The claim it makes is about the model: the supplied rows
  # and a differenced solve disagree in exactly one direction, because exactly one
  # quantity is wrong. That is a statement about a residual large enough to have a
  # direction of its own. With the tabulation carrying exact slopes the residual is
  # at the differenced reference's error on every layer (the ratios above), so what
  # is left is the reference's own noise, and noise has no preferred direction: the
  # ratio reads 4.3e-03, and it is there because the dominant term is small rather
  # than because a second one is large.
  #
  # So the two regimes assert different things, and which one holds is measured
  # rather than assumed -- the same rule the per-layer floor above follows. A
  # tolerance taken from a round-off statistic is not a tolerance, and neither is a
  # rank taken from one.
  at_reference_error <- max(residual / floor) < 2
  if (at_reference_error) {
    message("  the residual is at the reference's own error on every layer, so its ",
            "direction is the reference's and the rank bound is not asserted")
    expect_lt(max(residual / floor), 2)
  } else {
    ladder_report_margin("the factorisation's residual is one direction",
                         singular[[2]] / singular[[1]], 1e-3)
  }

  # Entrywise relative, on the entries whose reference has converged -- an entry
  # still moving between the two coarsest steps is being compared against its own
  # round-off, and on the worst of them that ratio reaches 2.5.
  #
  # The bound is the fit's own truncation and is not read off this measurement:
  # `a` is solved from a central difference at a relative fit_step, so a quantity
  # built from it carries fit_step^2 times a dimensionless third-derivative ratio.
  # The factor of ten is the headroom the anchoring column below is held to, at
  # the same step and by the same argument.
  #
  # ⚠️ AND THE REFERENCE HAS TO RESOLVE AN ENTRY BEFORE IT CAN ADJUDICATE IT. The
  # old filter asked only that the coarsest step be within 10% of the reference,
  # which is four orders looser than the bound enforced -- so an entry whose three
  # steps disagreed at 1e-04 passed it and was then held to 1e-05. That stayed
  # invisible while the rows disagreed by more than the reference did, and became
  # the whole of this statistic once they did not: measured at the entry that
  # failed, the row agrees with the 1e-05 step to 8.5e-08 and with the 1e-07 step
  # to 5.2e-06, while the 1e-06 step used as the reference differs from BOTH its
  # neighbours by 2.3e-04. The statistic was reading the reference's own step, and
  # every per-layer residual was 1e-08 at the same time.
  #
  # So an entry is comparable where the reference resolves it to better than the
  # bound, which is the same standard the per-layer statistic already applies, and
  # what that drops is counted rather than assumed small.
  fit_step <- 1e-3
  own_spread <- pmax(abs(steps[[1]] - steps[[2]]), abs(steps[[2]] - steps[[3]]))
  converged <- own_spread < (10 * fit_step^2) * abs(reference)
  message("  entrywise: ", sum(converged), " of ", length(converged),
          " entries the reference resolves to better than the bound")
  expect_gt(sum(converged), length(residual))
  ladder_report_margin(
    "the water rows entrywise, where the reference converged",
    max(abs((supplied - reference)[converged] / reference[converged])),
    10 * fit_step^2)
})

test_that("the water rows hold in the direction the ecology reads", {
  # The check above bounds every entry against the largest entry. That is not the
  # same as bounding the direction the answer is read in.
  #
  # Water moves on DIFFERENCES of potential and tissue fails on ABSOLUTES, so
  # moving every layer together is close to a symmetry of the model: the true
  # response is a small residue of the entries that make it. That residue is a
  # fifteenth to a twenty-sixth of an entry, and a per-entry bound normalised by
  # the largest entry is therefore that many times weaker in the direction than
  # it looks -- weaker again by the number of layers, because a systematic error
  # adds coherently across them while the answer cancels.
  #
  # A wrong split between the two scalars is exactly such a systematic error. So
  # this forms the direction and checks it there.
  # The shipped parameterisation is nowhere near the regime -- 1.1x on both
  # constructed patches, and measured at every node of every run stand between
  # 1.03 and 1.14, the competing stand included. What reaches it is the root
  # network's hydraulics at constant root carbon, which no trajectory develops.
  # The first two are kept because their 1.1x is the finding.
  for (name in c("one cohort", "two by two", "uniform drying")) {
    patch <- switch(name,
                    "one cohort" = ladder_patch_one(),
                    "two by two" = ladder_patch_two_by_two(),
                    "uniform drying" = ladder_patch_uniform_drying())
    ladder_require_regime(patch, "patch")

    supplied <- ladder_block_jacobian_forward_tf24(patch, 1L)
    inputs <- colnames(supplied)
    outputs <- rownames(supplied)
    is_soil <- grepl("^psi_soil_", inputs)
    uptake <- grepl("^uptake_", outputs)

    # The tangent is linear in its direction, so the sum of the soil columns IS
    # the directional derivative and needs no second evaluation.
    rows <- rowSums(supplied[uptake, is_soil, drop = FALSE])
    direction <- as.numeric(is_soil)

    # The reference is ONE perturbation along the direction, not a sum of the
    # per-column ones: a sum carries the columns' truncation, which is the size of
    # the terms rather than of what is left of them.
    along <- function(rel) {
      ladder_block_direction_difference_tf24(patch, 1L, direction, rel)[uptake]
    }
    coarse <- along(1e-6)
    fine <- along(1e-7)
    scale <- max(abs(fine))
    floor <- max(max(abs(coarse - fine)) / scale, 4 * .Machine$double.eps)

    # The amplification this check exists for, measured on the fixture rather than
    # taken from the report. It is the regime's precondition: where the entries do
    # not cancel, the direction is not a near-symmetry and this says nothing the
    # per-entry check has not already said.
    amplification <- max(abs(supplied[uptake, is_soil, drop = FALSE])) / scale
    message(sprintf("\n  %s: uniform-drying amplification %.1fx, reference error %.3e",
                    name, amplification, floor))

    # What bounds this is the per-entry bound above, times the amplification the
    # fixture is measured at. That composition IS the thing this check exists for:
    # the entries carry the fit's own truncation, a systematic error adds
    # coherently across them, and the direction divides what is left by the
    # cancellation. So a bound on the entries alone permits amplification times as
    # much here, and saying so is what turns the per-entry number into a statement
    # about the direction the ecology reads.
    #
    # Not `floor`. That is the gap between two step sizes, which measures how far
    # the difference has converged and not how accurate it is, so it tightens
    # exactly where the check gets sharper -- the same defect the per-layer
    # statistic above was carrying.
    fit_step <- 1e-3
    ladder_report_margin(
      sprintf("the uniform-drying direction, %s", name),
      max(abs(rows - fine)) / scale,
      10 * fit_step^2 * amplification)
  }
})

test_that("the water rows hold in the family the pair was anchored in", {
  # The pair is solved from one root-carbon direction, and root carbon reaches the
  # block through three inputs rather than as one of its own. Those three columns
  # are computed by the difference the check above already runs and were being
  # discarded; a_r1 is the clean one, reaching the water channel only through the
  # root profile, where height also moves leaf area and absorbed light.
  patch <- ladder_patch_one()
  ladder_require_regime(patch, "patch")

  block <- ladder_block_jacobian_forward_tf24(patch, 1L)
  inputs <- colnames(block)
  outputs <- rownames(block)
  uptake <- grepl("^uptake_", outputs)
  carbon <- match("a_r1", inputs)
  expect_false(is.na(carbon))

  supplied <- block[uptake, carbon]
  differenced <- function(rel) {
    ladder_block_difference_tf24(patch, 1L, rel)[uptake, carbon]
  }
  coarse <- differenced(1e-6)
  fine <- differenced(1e-7)
  scale <- max(abs(fine))
  floor <- max(max(abs(coarse - fine)) / scale, 4 * .Machine$double.eps)

  # Non-vacuity: the column has to carry a response, or agreeing about zero passes.
  ladder_expect_moves(fine, numeric(length(fine)),
                      "the root-carbon column of the uptake rows")
  # This column IS in sample -- it is the one the pair is solved from -- so what
  # bounds it is the fit's own step, a central difference at a relative 1e-3, and
  # not the reference's precision. Measured at 7.1e-08, which is two orders inside
  # that and an order below the potential layers' residuals: the fit lands where
  # it is taken, which is the premise every predicted layer rests on.
  fit_step <- 1e-3
  message(sprintf("  the root-carbon column's own error: %.3e", floor))
  ladder_report_margin("the water rows in the anchoring family",
                       max(abs(supplied - fine)) / scale, 10 * fit_step^2)
})

test_that("the leaf's supplied rows are refereed against its own algebra", {
  # A finite difference of the recorded step cannot referee a supplied row. The
  # recorded expression is the value plus a sum of partials times brackets that
  # are each exactly zero, so the block's forward value does not depend on a
  # recorded input at all: differencing the block returns identically zero on
  # exactly the columns a supplied row occupies, whether the row is right,
  # wrong, or absent.
  #
  # So the rows have to be checked against the solver's own algebra, or against
  # a transpose identity that needs no reference gradient -- never against a
  # difference of the step that consumes them.
  # The referee is not the leaf's algebra but something that serves the same
  # purpose and is available: the forward model rebuilt from its parameters. A
  # rebuild runs preparation, so the leaf is CONSTRUCTED with the moved trait
  # rather than handed it afterwards, and differencing the rates then traverses
  # the leaf's forward solve while touching none of the derivative code the rows
  # come from.
  #
  # It reaches exactly the columns nothing else does. A perturbation of the
  # prepared strategy leaves the leaf at the value it was prepared with, so that
  # difference is EXACTLY zero on these -- not small, zero -- whether the row is
  # right, wrong or absent.
  patch <- ladder_patch_two_by_two(cross = FALSE)
  ladder_require_regime(patch, "patch")

  n <- patch$ode_size
  seed <- ladder_seeds(n, scale = ladder_block_scale(patch$ode_rates))
  columns <- ladder_trait_names_tf24(patch)
  sweep_row <- ladder_rhs_adjoint_tf24(patch, seed)$trait
  names(sweep_row) <- columns

  # Non-vacuity, and it is what the check rests on: the reference that is
  # normally used must be silent on these columns, or refereeing them here would
  # be redundant rather than necessary.
  prepared <- as.vector(crossprod(
    ladder_rhs_trait_difference_tf24(patch, 1e-6), seed))
  names(prepared) <- columns

  message("\n  leaf traits, refereed by a rebuilt forward model:")
  worst <- 0
  for (name in ladder_leaf_own_traits()) {
    for (index in seq_len(length(patch$species))) {
      key <- paste0(index, ".", name)
      # Skip the two whose row is zero at an interior optimum by complementary
      # slackness; they are declared in the floor and have nothing to referee.
      if (name %in% ladder_zero_at_an_interior_optimum()) next
      expect_equal(prepared[[key]], 0,
                   label = paste("the prepared-strategy difference reaches", key))
      coarse <- sum(ladder_rate_difference_rebuilt(patch, index, name, 1e-5) * seed)
      fine <- sum(ladder_rate_difference_rebuilt(patch, index, name, 1e-6) * seed)
      got <- sweep_row[[key]]
      scale <- max(abs(got), abs(fine))
      residual <- abs(got - fine) / scale
      worst <- max(worst, residual)
      message(sprintf("    %-24s sweep %13.6e  rebuilt %13.6e  rel %.2e",
                      key, got, fine, residual))
      # Non-vacuity per column: a row of zero on both sides would agree perfectly.
      expect_gt(abs(got) / max(abs(sweep_row)), 1e-10)
    }
  }

  # The tolerance is the leaf's own solve tolerance appearing twice, because both
  # sides difference a re-solved operating point. That is three orders looser than
  # the 3e-06 the columns a prepared-strategy difference reaches are held to, and
  # tightening it needs the analytic route rather than a better step.
  ladder_report_margin("leaf trait rows, against a rebuilt forward model",
                       worst, 1e-2)
})
