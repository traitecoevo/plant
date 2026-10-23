# One cohort, and the Jacobian formed entirely.
#
# Use the strongest referee the size permits, and the size permits an exhaustive
# one exactly once -- here, at the bottom. Nothing above this is small enough to
# form whole, so nothing above it is ever checked this well again.
#
# Do not contract. A dot-product identity yields one number per run: an error in
# one cell enters as a seed component times that cell times an input component,
# so it hides behind a small seed, it cancels against other cells, and when it
# does fail it localises to nothing. At this size the whole object can be formed,
# so forming it is strictly better -- a failure names an output row and an input
# column rather than reporting that a sum disagreed.
#
# The unit is the adjoint of one right-hand-side evaluation at one stage. No
# solver, no trajectory.

test_that("the recorded block has the shape the design is costed on", {
  patch <- ladder_patch_one()
  block <- ladder_block_jacobian_forward_tf24(patch, 1L)
  inputs <- colnames(block)
  outputs <- rownames(block)

  # Six own states, the field's knot values and slopes, the canopy top its knots
  # are fractions of, the soil potentials and the traits in; six rates, the
  # density rate and one draw per layer out. Both counts move with the
  # configuration and neither is a property of the model; the shape is. Every
  # count is read off the names for that reason: a parameter added to the
  # strategy moves them together, where a restated total fails here and says
  # nothing about the shape.
  #
  # ⚠️ THE CANOPY TOP IS AN INPUT AND WAS NOT, which is what let the field's grid
  # sever a channel nothing could see. The knots sit at u_k * height_max, so a
  # block that reads values and slopes and not the top reads every knot and not
  # the grid they are on.
  n_knot <- sum(grepl("^light_value_", inputs))
  n_layer <- sum(grepl("^psi_soil_", inputs))
  n_par <- length(inputs) - 6L - 2L * n_knot - 1L - n_layer
  expect_equal(sum(grepl("^light_slope_", inputs)), n_knot)
  expect_equal(sum(inputs == "light_height_max"), 1L)
  expect_gt(n_par, 0L)
  expect_equal(length(outputs), 6L + 1L + n_layer)
  expect_equal(head(inputs, 6L),
               c("height", "mortality", "fecundity", "area_heartwood",
                 "mass_heartwood", "storage"))
  expect_equal(tail(outputs, 5L), paste0("uptake_", 1:5))
  message(sprintf(
      "\n  block shape: %d inputs (%d knots, %d layers, %d traits) by %d outputs",
      length(inputs), n_knot, n_layer, n_par, length(outputs)))
})

test_that("the block's forward and reverse Jacobians agree entry by entry", {
  patch <- ladder_patch_one()
  ladder_require_regime(patch, "patch")

  forward <- ladder_block_jacobian_forward_tf24(patch, 1L)
  reverse <- ladder_block_jacobian_reverse_tf24(patch, 1L)
  inputs <- colnames(forward)
  outputs <- rownames(forward)

  # Non-vacuity: a Jacobian of a block that returned nothing is not evidence.
  expect_true(all(is.finite(forward)))
  expect_true(all(is.finite(reverse)))
  expect_gt(sum(forward != 0), length(outputs))

  # Each entry is measured against its own output row's scale, not against
  # itself. A cell whose true value is a hundred orders below its row carries no
  # information about the row, and a per-cell relative measure would make the
  # tolerance a statement about that cell's round-off rather than about the
  # Jacobian. The scale is the larger of the two rows, which is the one statistic
  # ladder_matrix_residual carries, so the corruption tests measure this comparison
  # rather than a second one written to resemble it.
  row_scale <- pmax(apply(abs(forward), 1, max), apply(abs(reverse), 1, max),
                    .Machine$double.xmin)
  residual <- abs(forward - reverse) / row_scale
  worst <- which(residual == max(residual), arr.ind = TRUE)[1, ]
  message(sprintf("  worst cell: %s / %s  forward %.6e  reverse %.6e",
                  outputs[[worst[[1]]]], inputs[[worst[[2]]]],
                  forward[worst[[1]], worst[[2]]],
                  reverse[worst[[1]], worst[[2]]]))
  expect_equal(max(residual), ladder_matrix_residual(reverse, forward))

  ladder_report_margin("block Jacobian, forward against reverse",
                       max(residual), 10 * ladder_forward_floor(patch))
})

test_that("the block's structure is what every cost argument assumes", {
  # The matrix formed above is the five-way classification, measured rather than
  # asserted. Asserting it directly costs a rank computation on a matrix that
  # already exists, and this is the one place the structure is checked rather
  # than read.
  patch <- ladder_patch_one()
  j <- ladder_block_jacobian_forward_tf24(patch, 1L)
  inputs <- colnames(j)
  outputs <- rownames(j)

  is_light <- grepl("^light_", inputs)
  is_value <- grepl("^light_value_", inputs)
  is_slope <- grepl("^light_slope_", inputs)
  is_soil <- grepl("^psi_soil_", inputs)
  uptake <- grepl("^uptake_", outputs)

  numerical_rank <- function(m, tol = 1e-9) {
    if (all(m == 0)) return(0L)
    s <- svd(m)$d
    sum(s > tol * max(s))
  }

  # The field reaches the physiology through one scalar, so the whole knot block
  # is rank one on this coordinate. On the height coordinate the density rate
  # would evaluate the physiology a second time at a displaced height and make
  # its own field query, and the block would have rank at least two.
  light_block <- j[, is_light, drop = FALSE]
  message(sprintf("  light block rank %d over %d columns",
                  numerical_rank(light_block), sum(is_light)))
  expect_equal(numerical_rank(light_block), 1L)

  # Per-layer uptake reads its own layer's potential and the collar, so the block
  # is a diagonal explicit part plus the argmax channel, and the argmax channel is
  # rank one because the collar is one number: M = D + u v'.
  #
  # The explicit part is NOT the observed diagonal. Zeroing the diagonal removes
  # D + diag(u_i v_i), which takes the rank-one term's own diagonal with it and
  # leaves a matrix of full rank -- so a rank test on that residue fails whether
  # or not the claim holds, and measuring it here gives five singular values of
  # comparable size rather than one.
  #
  # What is testable without knowing D is the OFF-DIAGONAL structure, and it is
  # testable as a prediction rather than a fit. For M_ij = u_i v_j off the
  # diagonal, one reference row and one reference column determine u and v, and
  # every remaining off-diagonal entry follows:
  #
  #   M_ij = M_ic * M_rj / M_rc,   i != j, i != c, j != r
  #
  # None of the predicted entries is used to build the prediction, so this is out
  # of sample, and the reference entry's magnitude is reported beside it so a
  # degenerate pair is refused rather than absorbed.
  soil_block <- j[uptake, is_soil, drop = FALSE]
  expect_equal(dim(soil_block), c(5L, 5L))
  n_layer <- nrow(soil_block)
  scale <- max(abs(soil_block))

  # The diagonal is the explicit Ohm's-law term and must dominate: a layer's own
  # potential drives its own flux. Without this the split is not identifiable and
  # the prediction below would be a statement about noise.
  expect_true(all(abs(diag(soil_block)) >
                    apply(abs(soil_block) - diag(diag(abs(soil_block))), 1, max)))

  r <- 1L; c <- 2L
  pivot <- soil_block[r, c]
  expect_gt(abs(pivot) / scale, 1e-6)
  predicted <- 0
  worst <- 0
  for (i in seq_len(n_layer)) {
    for (k in seq_len(n_layer)) {
      if (i == k || i == c || k == r) next
      got <- soil_block[i, k]
      want <- soil_block[i, c] * soil_block[r, k] / pivot
      predicted <- predicted + 1L
      worst <- max(worst, abs(got - want) / max(abs(got), abs(want)))
    }
  }
  message(sprintf(
    "  uptake-by-potential: rank %d, off-diagonal rank one predicted over %d entries",
    numerical_rank(soil_block), predicted))
  expect_gt(predicted, 0L)
  ladder_report_margin("argmax channel, out-of-sample rank one",
                       worst, 1e4 * ladder_forward_floor(patch))

  # Density enters no cohort rate: it reaches the world only through the two
  # reductions. It is not among the block's inputs at all, which is the strongest
  # form of that claim.
  expect_false(any(grepl("log_density", inputs)))

  # The value and slope channels are separate inputs and both must be live. An
  # identity posed on values alone passes with the entire slope channel absent.
  ladder_expect_moves(j[, is_value], matrix(0, nrow(j), sum(is_value)),
                      "the field's value channel")
  ladder_expect_moves(j[, is_slope], matrix(0, nrow(j), sum(is_slope)),
                      "the field's slope channel")
})

test_that("every trait the block reads has a column, or is refused by name", {
  # An absent row is an exact zero that reads as an ecological finding. The block
  # is where a trait's row either exists or does not, so this is where the
  # classification is cheapest to make.
  patch <- ladder_patch_one()
  j <- ladder_block_jacobian_forward_tf24(patch, 1L)
  inputs <- colnames(j)
  traits <- inputs[!grepl("^light_|^psi_soil_", inputs)][-(1:6)]
  columns <- j[, match(traits, inputs), drop = FALSE]
  peak <- apply(abs(columns), 2, max)
  names(peak) <- traits

  refused <- ladder_refused_by_name()
  # This level has its own two causes and not the census's. The parameters the
  # census cannot see for want of a metric reading offspring are live HERE,
  # because the offspring-production rate is one of the block's own outputs, so
  # excusing them here would excuse an absent row; and the parameters that reach
  # the census only through the introduction boundary are correctly silent here,
  # because the block is one individual's physiology and no introduction is in it.
  declared <- c(ladder_zero_at_an_interior_optimum(),
                ladder_zero_outside_the_cohort_block(),
                ladder_zero_outside_every_rate())
  unaccounted <- names(peak)[peak == 0 & !(names(peak) %in% c(refused, declared))]
  if (length(unaccounted) > 0) {
    message("\n  traits with no row in the block: ",
            paste(unaccounted, collapse = " "))
  }
  expect_length(unaccounted, 0L)

  # And the pair the census cannot see is live here, which is what makes that a
  # statement about the metric set rather than about the model.
  for (name in ladder_zero_outside_the_metric_support()) {
    expect_gt(peak[[name]], 0)
  }
})

# The reduction transposes, where the reference shares no code.
#
# The checks below run on the uncrossed fixture, and that is a concession rather
# than a choice: the crossed one refuses, and the test directly above owns that
# refusal. Once the transposes test the abscissa rather than the height, these
# should move to the crossed fixture, which is where the design says every run
# belongs.

test_that("a crossed stand transposes, because the forward model runs it", {
  # Reserve-gated growth lets a younger cohort overtake an older one, so on the
  # birth-date coordinate crossing is more common, not less. The quantity that
  # must stay monotone here is the abscissa, and the abscissa is not height -- so
  # a guard written as a test of height ordering refuses exactly the stands the
  # forward model handles correctly.
  #
  # Both failures are silent and they are silent in opposite directions: a
  # reduction that integrates a crossed grid without guarding has neighbouring
  # trapezia cancelling instead of accumulating, and a guard on the wrong
  # quantity stops on a stand that is fine.
  crossed <- ladder_patch_two_by_two(cross = TRUE)
  ladder_require_regime(crossed, "patch")
  nodes <- ladder_nodes(crossed)
  expect_gt(ladder_crossing_count(nodes), 0L)

  # The forward model runs it, which is the half of the claim that has to hold
  # before the other half is a defect rather than a preference.
  rates <- crossed$ode_rates
  expect_true(all(is.finite(rates)))

  seed <- ladder_seeds(crossed$ode_size,
                       scale = ladder_block_scale(rates))
  expect_no_error(ladder_rhs_adjoint_tf24(crossed, seed))
})

test_that("the tangent carries the soil channel, and two references agree in it", {
  # The environment holds its integrated soil state at the scalar the model
  # carries, so a tangent seeded anywhere propagates through the water balance and
  # back through the retention curve. Every soil column and every soil rate row of
  # its Jacobian is non-zero, and that is asserted rather than left to a
  # residual, because a reference silent about a channel and a transpose wrong in
  # it agree perfectly.
  patch <- ladder_patch_two_by_two(cross = FALSE)
  ladder_require_regime(patch, "patch")

  n <- patch$ode_size
  n_env <- patch$environment$ode_size
  is_soil <- seq_len(n) > (n - n_env)
  jacobian <- ladder_rhs_state_jacobian_forward_tf24(patch)
  expect_equal(dim(jacobian), c(n, n))
  expect_true(all(is.finite(jacobian)))

  expect_gt(sum(jacobian[, is_soil] != 0), n_env)
  expect_gt(sum(jacobian[is_soil, ] != 0), n_env)

  # The two references are independent -- one exact and tape-free, the other a
  # plain-double difference that re-runs the forward balance -- so their agreement
  # in the channel that was unrefereed is what licenses using either in it.
  coarse <- ladder_rhs_state_difference(patch, 1e-5)
  fine <- ladder_rhs_state_difference(patch, 1e-6)
  floor <- ladder_difference_floor(coarse, fine)
  ladder_report_margin("tangent against a difference, whole Jacobian",
                       max(abs(jacobian - fine)) / max(abs(fine)),
                       1e3 * floor)

  # The transpose against the tangent, seeded on every rate. Off the soil rates it
  # is round-off. On them it is the leaf's own supplied-row gap: both reductions
  # now carry the boundary node's row, so what is left is the difference between a
  # recorded uptake row and one taken by re-solving the collar.
  rates <- patch$ode_rates
  seed <- ladder_seeds(n, scale = ladder_block_scale(rates))
  cohort_seed <- seed
  cohort_seed[is_soil] <- 0
  expected <- as.vector(crossprod(jacobian[!is_soil, , drop = FALSE],
                                  cohort_seed[!is_soil]))
  observed <- ladder_rhs_adjoint_tf24(patch, cohort_seed)$state
  ladder_report_margin("right-hand-side transpose, off the soil rates",
                       max(abs(expected - observed)) / max(abs(expected)),
                       1e3 * ladder_forward_floor(patch))

  soil_seed <- seed
  soil_seed[!is_soil] <- 0
  soil_expected <- as.vector(crossprod(jacobian[is_soil, , drop = FALSE],
                                       seed[is_soil]))
  soil_observed <- ladder_rhs_adjoint_tf24(patch, soil_seed)$state
  ladder_report_margin("the transpose on the soil rates",
                       max(abs(soil_expected - soil_observed)) /
                         max(abs(soil_expected)),
                       3 * ladder_soil_row_agreement())
})

test_that("one right-hand-side evaluation transposes to a difference of itself", {
  # The referee for the channel the tangent cannot reach. A plain-double central
  # difference of the same right-hand side traverses the forward soil balance and
  # the retention curve, and the transpose under test is not on its path -- which
  # is what path disjointness asks for, and it is not the same requirement as
  # sharing no code.
  #
  # A difference is admissible HERE and not at the block, and the distinction is
  # the supplied row: the leaf's rows are recorded, so a difference of the step
  # that consumes them is identically zero on those columns whether they are
  # right, wrong or absent. The soil balance records nothing, so differencing it
  # measures what it computes. This is where d(psi)/d(theta) and the drainage
  # cascade's own derivative get a referee for the first time.
  patch <- ladder_patch_two_by_two(cross = FALSE)
  ladder_require_regime(patch, "patch")

  n <- patch$ode_size
  rates <- patch$ode_rates
  seed <- ladder_seeds(n, scale = ladder_block_scale(rates))
  observed <- ladder_rhs_adjoint_tf24(patch, seed)$state

  # Two steps, so the tolerance is the reference's own uncertainty rather than a
  # literal: a central difference is truncation-dominated above its minimum and
  # round-off-dominated below, so the gap between two decades bounds the better
  # one's error.
  coarse <- as.vector(crossprod(ladder_rhs_state_difference(patch, 1e-5), seed))
  fine <- as.vector(crossprod(ladder_rhs_state_difference(patch, 1e-6), seed))
  floor <- ladder_difference_floor(coarse, fine)
  message(sprintf("  the difference's own error at this step: %.3e", floor))

  scale <- max(abs(fine))
  residual <- max(abs(fine - observed)) / scale
  worst <- which.max(abs(fine - observed))
  message(sprintf("  worst state entry %d: difference %.6e  adjoint %.6e",
                  worst, fine[[worst]], observed[[worst]]))

  # Non-vacuity: the soil block has to carry adjoint traffic, or this check
  # passes with the whole water channel missing from both sides.
  n_env <- patch$environment$ode_size
  is_soil <- seq_len(n) > (n - n_env)
  ladder_expect_moves(observed[is_soil], rep(0, sum(is_soil)),
                      "the soil block's state adjoint")

  ladder_report_margin("right-hand-side transpose, against a difference",
                       residual, 10 * floor)
})

test_that("the state Jacobian holds at the states a trajectory reached", {
  # Every check above forms its Jacobian on a patch whose state was written by
  # hand and conditioned into the declared regime. A trajectory's state is
  # wherever the run went, and "correct at one state" and "correct at every state
  # the trajectory visits" are different claims. This is the second.
  #
  # Forward against reverse, and not against a difference. A difference of the
  # right-hand side re-solves the leaf while both differentiated paths consume its
  # supplied rows, so a difference here prices the factorisation's own fit rather
  # than the transpose: measured flat in the step and landing on the uptake
  # accumulator, which is the documented soil-row gap. Two exact paths cannot hide a
  # supplied row behind each other, so the whole matrix is comparable at the
  # arithmetic floor.
  stand <- ladder_stand_trajectory(lifetime = 2)
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  # The first recorded state is the one no step reached, and on this coordinate it
  # holds the environment and no cohort. From the introduction onward the width is
  # fixed, which is what lets one patch load every state.
  expect_equal(length(unique(widths[-1])), 1L)

  patch <- ladder_as_patch(stand)
  names_state <- ladder_rate_names(patch)
  n <- patch$ode_size
  at <- unique(round(seq(2, length(trajectory), length.out = 6)))
  heights <- numeric(0)
  worst_overall <- 0
  worst_label <- ""

  for (k in at) {
    patch$set_ode_state(trajectory[[k]]$state, trajectory[[k]]$time)
    # A state load leaves the inflow condition at its FIRST evaluation, taken with
    # every species' boundary interval left off; a run carries the second, taken in
    # the field rebuilt to include it. Asking for the rates is what advances it from
    # R, where set_recorded_state is not reachable.
    #
    # Without this the two paths linearise at different boundary nodes and disagree
    # by their own whole magnitude -- 1.0 at every loaded state, against 3e-16 with
    # it. That is an instrument reading, and it is the reason this line exists
    # rather than a tolerance being widened to cover it.
    invisible(patch$ode_rates)

    forward <- ladder_rhs_state_jacobian_forward_tf24(patch)
    # One unit output adjoint per rate gives one row of the transpose, so the whole
    # object is formed rather than contracted: a contraction against one direction
    # lets a wrong cell hide behind a small seed component and localises to nothing.
    reverse <- t(vapply(seq_len(n), function(i) {
      ladder_rhs_adjoint_tf24(patch, replace(numeric(n), i, 1))$state
    }, numeric(n)))

    # Scaled by the larger of the two rows. Where a row is near zero on one path and
    # not on the other, dividing by the small one reports a ratio of two round-offs
    # rather than a disagreement.
    row_scale <- pmax(apply(abs(forward), 1, max), apply(abs(reverse), 1, max),
                      .Machine$double.xmin)
    resid <- abs(forward - reverse) / row_scale
    cell <- which(resid == max(resid), arr.ind = TRUE)[1, ]
    height <- patch$species[[1]]$nodes[[1]]$height
    heights <- c(heights, height)
    message(sprintf("    state %3d  h %6.3f  worst %.2e at %s by %s", k, height,
                    max(resid), names_state[[cell[[1]]]],
                    names_state[[cell[[2]]]]))
    if (max(resid) > worst_overall) {
      worst_overall <- max(resid)
      worst_label <- sprintf("state %d, h %.3f, %s by %s", k, height,
                             names_state[[cell[[1]]]], names_state[[cell[[2]]]])
    }
  }

  # Non-vacuity, twice over: the sampled states have to be different states, and
  # the matrix has to have entries.
  expect_gt(diff(range(heights)), 1)
  expect_gt(sum(forward != 0), n)
  message(sprintf("  worst over the trajectory: %.2e at %s", worst_overall,
                  worst_label))
  # A hundred times the fixture's floor rather than ten. These states are loaded
  # rather than constructed and the field is rebuilt from a recorded vector, so a
  # few extra machine units are the cost of the load; measured worst is 7.6e-15
  # against a floor of 8.9e-16.
  ladder_report_margin("state Jacobian along a trajectory, forward against reverse",
                       worst_overall, 100 * ladder_forward_floor(patch))
})

test_that("the trait rows that arise inside a reduction arrive", {
  # Four parameters reach a census through the field build as well as through a
  # plant: the extinction coefficient, the crown shape, and the two allometric
  # constants that set leaf area at fixed height. A transpose written through a
  # structure holding only size and density slots has no route to a trait
  # accumulator however correct its arithmetic, and those rows then read exactly
  # zero.
  patch <- ladder_patch_two_by_two(cross = FALSE)

  n <- patch$ode_size
  n_env <- patch$environment$ode_size
  is_soil <- seq_len(n) > (n - n_env)
  columns <- ladder_trait_names_tf24(patch)
  rates <- patch$ode_rates
  seed <- ladder_seeds(n, scale = ladder_block_scale(rates))

  # Seeded on the cohort rates only. The soil rates now carry trait rows in the
  # reference too, but their transpose is short of the boundary node's own draw,
  # so including them here would report that gap as a reduction defect.
  cohort_seed <- seed
  cohort_seed[is_soil] <- 0
  jacobian <- ladder_rhs_trait_jacobian_forward_tf24(patch)
  expect_gt(sum(jacobian[is_soil, ] != 0), n_env)
  expected <- as.vector(crossprod(jacobian[!is_soil, , drop = FALSE],
                                  cohort_seed[!is_soil]))
  observed <- ladder_rhs_adjoint_tf24(patch, cohort_seed)$trait
  names(expected) <- names(observed) <- columns

  # Non-vacuity first: the reduction's own contribution has to be non-zero in
  # the reference, or its absence in the sweep proves nothing. Both species'
  # columns are taken, because a per-species row summed over every cohort of
  # every species is the collapse this fixture exists to catch. The crown shape
  # is the fourth reduction-borne parameter and has no column at all, so it drops
  # out here rather than being checked.
  reduction_borne <- c("k_I", "eta", "a_l1", "a_l2")
  present <- names(expected)[ladder_bare_traits(names(expected)) %in%
                               reduction_borne]
  expect_length(present, 6L)
  expect_true(all(abs(expected[present]) > 0),
              label = "the tangent's reduction-borne rows are themselves zero")

  scale <- max(abs(expected))
  gap <- names(expected)[abs(expected - observed) > 1e-6 * scale]
  if (length(gap) > 0) {
    message("\n  trait rows that disagree: ", paste(unique(gap), collapse = " "))
  }
  ladder_report_margin("trait transpose, off the soil rows",
                       max(abs(expected - observed)) / scale,
                       1e3 * ladder_forward_floor(patch))
})

test_that("the trait rows a difference can reach agree with one", {
  # The whole contraction, soil rows included, refereed against a plain-double
  # difference of the same right-hand side in each trait. This is the only
  # reference that carries the water channel, so it is the only one that referees
  # a trait whose route runs through it.
  #
  # Eighteen columns are outside its reach and the set is declared rather than
  # discovered: the leaf holds its own copy of nine traits per species, taken when
  # the strategy was prepared, and a rate evaluation does not push them back in.
  # So the difference is EXACTLY zero there whether the sweep's row is right,
  # wrong or absent -- the same shape as a recorded row one level down. Their
  # referee is the forward model rebuilt from its parameters, which
  # `test-gradient-ladder-factorisation.R` applies to exactly this set of traits.
  patch <- ladder_patch_two_by_two(cross = FALSE)

  n <- patch$ode_size
  columns <- ladder_trait_names_tf24(patch)
  seed <- ladder_seeds(n, scale = ladder_block_scale(patch$ode_rates))
  observed <- ladder_rhs_adjoint_tf24(patch, seed)$trait

  coarse <- as.vector(crossprod(
    ladder_rhs_trait_difference_tf24(patch, 1e-5), seed))
  fine <- as.vector(crossprod(
    ladder_rhs_trait_difference_tf24(patch, 1e-6), seed))

  # Which columns the difference cannot reach, asserted to be exactly the
  # declared set. A column that leaves this set has become refereeable, and one
  # that joins it has stopped being so; either is a finding.
  scale_all <- max(abs(fine))
  unreachable <- columns[fine == 0 & abs(observed) > 1e-8 * scale_all]
  message("\n  columns no difference of the rates can reach: ",
          paste(unique(ladder_bare_traits(unreachable)), collapse = " "))
  expect_setequal(ladder_bare_traits(unreachable), ladder_leaf_own_traits())

  keep <- setdiff(seq_along(columns), match(unreachable, columns))

  # A second set this reference cannot referee, on a different mechanism from the
  # declared one above: the eight that reach birth size. The difference moves a
  # prepared strategy in place, and on the double path the seed's height is the
  # value preparation already resolved -- so it holds birth size still while the
  # adjoint derives it from the residual.
  #
  # The model supplies its own control for that claim, which is why this is a
  # split rather than a loosened tolerance. Fecundity reads (omega + a_f3) as a
  # sum, so those two share a rate route exactly, and a_f3 does not enter the
  # seed's residual. Their differences agree to every digit printed; omega's
  # disagreement is an order larger and a_f3's is where every other column sits.
  # Nothing but the seed-height route separates them.
  birth <- ladder_bare_traits(columns) %in% ladder_birth_size_parameters()
  refereeable <- setdiff(keep, which(birth))
  floor <- ladder_difference_floor(coarse[refereeable], fine[refereeable])
  message(sprintf("  the difference's own error at this step: %.3e", floor))

  scale <- max(abs(fine[keep]))
  resid <- abs(fine - observed) / scale
  worst <- keep[[which.max(resid[keep])]]
  message(sprintf("  worst refereeable column %s: difference %.6e  adjoint %.6e",
                  columns[[worst]], fine[[worst]], observed[[worst]]))
  message(sprintf("  birth-size columns %.3e, every other column %.3e",
                  max(resid[intersect(keep, which(birth))]),
                  max(resid[refereeable])))
  # Earning the split: the largest disagreement has to be one of the eight, or it
  # is a defect that has nothing to do with birth size and the exclusion hides it.
  expect_true(birth[[worst]])
  ladder_report_margin(
    sprintf("trait transpose against a difference, %d columns",
            length(refereeable)),
    max(resid[refereeable]), 10 * floor)
})

test_that("the recording does not grow with the stand", {
  # Peak is one cohort's rates at one stage, so it is flat in the cohort count
  # and in the number of differentiation targets: seeding more traits adds
  # registered slots, not recorded operations. Either climbing with the cohort
  # count is the tape leaking.
  one <- ladder_patch_one()
  four <- ladder_patch_two_by_two()

  seed_one <- rep(1, one$ode_size)
  seed_four <- rep(1, four$ode_size)
  a <- ladder_rhs_adjoint_tf24(one, seed_one)
  b <- ladder_rhs_adjoint_tf24(four, seed_four)
  per_unit <- (b$block_recording_size - a$block_recording_size) / 3
  message(sprintf(
    "\n  stage recording: %.0f at 1 cohort, %.0f at 4, %.0f per unit",
    a$block_recording_size, b$block_recording_size, per_unit))

  # The stage is recorded whole, so peak holds every unit at once and the size is
  # linear in the unit count. That is the commitment traded away: a recording flat
  # in the cohort count was paid for by writing every transpose between the state
  # and the rates by hand.
  #
  # What replaces a flatness test is amortisation. Four units cost less than four
  # one-unit recordings, because the field is read into the tape once for the
  # stage where a per-unit block registered the whole 130-wide field read again
  # for every unit. A recording that lost that would land at four times the one.
  expect_gt(per_unit, 0)
  expect_lt(b$block_recording_size, 4 * a$block_recording_size)
})
