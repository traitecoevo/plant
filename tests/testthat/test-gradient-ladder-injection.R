# Whether the checks bite, established by breaking what they watch.
#
# A suite that names an injected fault and never invokes it has switches on
# paper: every margin recorded elsewhere says how much room a check had, and none
# of them says the check would have noticed. Counting the
# injections is what separates a rung that has been climbed from one that has been
# written down.
#
# What these establish and what they do not. Each corruption is applied to one side
# of a comparison the suite already makes, so what is measured is the COMPARISON's
# sensitivity -- that a transpose broken in this named way would be rejected, and by
# what factor. None of them establishes that the C++ transpose is free of the
# defect; that is what the comparison itself does, on the uncorrupted object, in the
# rung that owns it. The two halves are worth nothing separately.
#
# Every injection here runs on one constructed patch and forms no trajectory, so the
# whole file is a second or two: a check whose sensitivity is expensive to establish
# does not get established.

test_that("the block Jacobian's comparison rejects each named corruption", {
  patch <- ladder_patch_one()
  forward <- ladder_block_jacobian_forward_tf24(patch, 1L)
  reverse <- ladder_block_jacobian_reverse_tf24(patch, 1L)
  inputs <- colnames(forward)
  outputs <- rownames(forward)
  bound <- 10 * ladder_forward_floor(patch)

  # The check as the rung makes it, so the injections below are measured against a
  # comparison that currently passes rather than against one that does not.
  ladder_report_margin("uncorrupted, forward against reverse",
                       ladder_matrix_residual(reverse, forward), bound)

  is_value <- grepl("^light_value_", inputs)
  is_slope <- grepl("^light_slope_", inputs)
  is_soil <- grepl("^psi_soil_", inputs)
  column <- function(name) which(inputs == name)

  # Each entry is one named defect, expressed as the arithmetic it would produce in
  # the transpose it breaks.
  corruptions <- list(
    # The slope channel's own contribution, which an identity posed on values alone
    # never touches: the pair is not a convenience, because the height transpose of
    # the slope channel carries a second-derivative term appearing nowhere in the
    # value transpose.
    "the slope channel dropped" =
      function(m) { m[, is_slope] <- 0; m },
    # A sign error survives any magnitude-only check, so it is worth its own row.
    "the slope channel's sign flipped" =
      function(m) { m[, is_slope] <- -m[, is_slope]; m },
    "the value channel dropped" =
      function(m) { m[, is_value] <- 0; m },
    # Half of the leaf-area-by-height product. The two halves are separately wrong
    # and jointly plausible, which is why dropping one is a listed fault.
    "one allometric constant's row dropped" =
      function(m) { m[, column("a_l1")] <- 0; m },
    # A reduction transpose written through a structure holding only size and
    # density slots has no route to a trait accumulator, however correct its
    # arithmetic. Its signature is a column of exact zeros.
    "a trait column routed nowhere" =
      function(m) { m[, column("k_I")] <- 0; m },
    # A trait treated as one input per cohort rather than one input read by every
    # cohort returns a fixed fraction of the right answer, with the correct sign and
    # no error raised. That is the shape to be sure of catching.
    "a trait column at a fixed fraction" =
      function(m) { m[, column("lma")] <- 0.5 * m[, column("lma")]; m },
    # The soil rows survive their own omission as a plausible bidiagonal transpose,
    # so their absence has to be caught by comparison rather than by inspection.
    "the soil potential rows dropped" =
      function(m) { m[, is_soil] <- 0; m },
    # The retention derivative applied twice rather than once. The row survives as
    # the potential's adjoint: dimensionally wrong, finite, and plausible.
    "the soil rows scaled twice" =
      function(m) { m[, is_soil] <- 2 * m[, is_soil]; m }
  )

  for (name in names(corruptions)) {
    broken <- corruptions[[name]](reverse)
    # A corruption that changed nothing would be detected at zero and read as a
    # check that does not work, when the fault is that the fixture has no such
    # channel to break.
    expect_false(identical(broken, reverse), label = paste("corruption is live:", name))
    ladder_report_detection(name, ladder_matrix_residual(broken, forward), bound)
  }
})

test_that("the structural assertions reject a structure that does not hold", {
  # The cost arguments in this corpus rest on four claims about the block's shape,
  # and the rung that owns them asserts each on the formed matrix. These are the
  # same claims broken on purpose: a rank test that passes on a matrix of any rank
  # is not a rank test.
  patch <- ladder_patch_one()
  j <- ladder_block_jacobian_forward_tf24(patch, 1L)
  inputs <- colnames(j)
  outputs <- rownames(j)

  numerical_rank <- function(m, tol = 1e-9) {
    if (all(m == 0)) return(0L)
    s <- svd(m)$d
    sum(s > tol * max(s))
  }
  is_light <- grepl("^light_", inputs)

  # Rank one over the whole knot block is what the birth-date coordinate buys. On
  # the height coordinate the density rate would evaluate the physiology a second
  # time at a displaced height and make its own field query, and the block would
  # have rank at least two -- so a second, independent column is exactly the defect
  # this claim excludes.
  light <- j[, is_light, drop = FALSE]
  expect_equal(numerical_rank(light), 1L)
  # A ramp rather than a rearrangement of an existing column: adding a multiple of
  # a column that is already in the range leaves the rank where it was, which is
  # what a rank-one matrix means.
  broken <- light
  ramp <- seq_len(nrow(broken)) / nrow(broken)
  broken[, 2] <- broken[, 2] + max(abs(light)) * ramp
  message(sprintf("  a second independent field column takes the rank from %d to %d",
                  numerical_rank(light), numerical_rank(broken)))
  expect_gt(numerical_rank(broken), 1L)

  # Density enters no cohort rate, and the strongest form of that claim is that it
  # is not among the block's inputs at all.
  expect_false(any(grepl("log_density", inputs)))
})
