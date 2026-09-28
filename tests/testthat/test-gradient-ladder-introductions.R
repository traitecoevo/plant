# Introductions, where the state changes dimension.
#
# The minimum fixture is two species and three introductions in the order
# species 1, species 2, species 1 -- four ranges, with the node stride
# exercised in both directions. Introducing into species 1 while species 2 exists
# is the arrangement a narrowing implemented as a truncation of the tail fails
# on.
#
# Each boundary channel is probed by a parameter reaching the census through that
# channel and no other, so a dropped channel gives an exact zero rather than a
# small error. That is the cleanest test shape in the ladder, and the boundary is
# where it is available.

test_that("the fixture actually widens, in both stride directions", {
  # Non-vacuity before anything else: a probe of introductions on a stand that
  # introduces nothing proves nothing.
  stand <- ladder_stand_introductions()
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  changes <- which(diff(widths) != 0)
  message(sprintf("\n  %d recorded steps, %d introductions, widths %s",
                  length(widths), length(changes),
                  paste(unique(widths), collapse = " -> ")))
  expect_gte(length(changes), 3L)

  # Species 1 is introduced into while species 2 already exists, which is what
  # makes a tail truncation wrong rather than merely narrow.
  nodes <- ladder_nodes(stand)
  first <- nodes[nodes$species == 1 & !nodes$boundary, ]
  second <- nodes[nodes$species == 2 & !nodes$boundary, ]
  expect_gt(max(first$birth_date), min(second$birth_date))
})

test_that("a reloaded state carries the boundary node the run carries", {
  # A stage evaluates the inflow condition twice, in two different fields, and
  # they are the same function at different arguments: the first is taken with
  # every species' boundary interval left off, the field is then rebuilt including
  # it, and the second is taken in that rebuilt field. The reductions were built
  # on the first; the water aggregation, an introduced node and a census all read
  # the second.
  #
  # Every rebuild the sweep does -- the census seed, the census's direct term, the
  # replay of an introduction, and the introduction's own transpose -- starts from
  # a recorded state, and a state load alone stops at the first evaluation. A sweep
  # that stops there linearises a boundary node the trajectory never carried, and
  # the only thing standing between that and a plausible wrong gradient is this
  # check: nothing about either number says which one a caller is holding.
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29, 0.94), c(0, 0.57))
  stand <- ladder_run(p)

  # What the run left behind is what the run introduced from and censused.
  left <- ladder_boundary_density(stand)
  patch <- ladder_as_patch(stand)
  both <- ladder_boundary_evaluations(patch, patch$ode_state, patch$ode_time)

  # Non-vacuity, and it is what the check rests on: if the two evaluations agreed
  # there would be no convention to get wrong and this check would pass on a
  # model that does not have the hazard.
  gap <- max(abs(both$in_uptake - both$in_field) / abs(both$in_field))
  message(sprintf("\n  the two evaluations differ by %.3e relative", gap))
  expect_gt(gap, 1e-6)

  # And the run carries the second, exactly. Tolerance is zero: this is a property
  # the implementation either has or does not.
  expect_identical(both$in_uptake, left)
  expect_false(isTRUE(all.equal(both$in_field, left, tolerance = 1e-12)))
})

test_that("the inflow condition's own parameters carry a row", {
  # The boundary node holds no ODE row, so its density adjoint has to be pushed
  # back through the condition that sets it: the adjoint at the boundary times the
  # boundary condition's own derivative. Accumulate it and never read it and the
  # accumulator is write-only -- every parameter whose only route to the census is
  # establishment then comes back short, finite and correctly signed.
  #
  # Two parameters reach the census through the inflow condition and nothing else,
  # which is what makes them the probe: the establishment constant, and the decay
  # of establishment in patch age. Neither appears in any cohort's rate equation.
  shared <- ladder_shared("introductions")
  stand <- shared$stand
  result <- shared$gradient
  columns <- colnames(result$gradient)

  for (name in c("1.a_d0", "2.a_d0", "1.recruitment_decay",
                 "2.recruitment_decay")) {
    reference <- ladder_trajectory_tangent(
      stand, ladder_trait_direction(columns, name))$tangent
    peak <- max(abs(result$gradient[, name]))
    expect_gt(peak, 0)
    ladder_report_margin(paste("column", name),
                         max(abs(reference - result$gradient[, name])) / peak,
                         ladder_trajectory_agreement())
  }
})

test_that("the boundary node's row carries the field its own cohorts build", {
  # A twin at an active scalar derives the auxiliaries a state determines as it is
  # built, which is before a caller has seeded anything. Leaf area at fixed height
  # then carries no row, and neither does the field the reduction builds from it,
  # so the two constants that set leaf area come back short by their whole field
  # channel while every parameter the reduction reads directly stays exact.
  #
  # The two groups are what localises it: k_I enters the reduction as its own
  # parameter, a_l1 and a_l2 only through a cohort's cached leaf area. Their sizes
  # are what makes it matter -- the seedling's own allometry and its chain through
  # the seed height cancel to about a tenth for this pair, so the omitted term is
  # larger than the residue and inverts the sign rather than shortening it.
  patch <- ladder_patch_one()
  for (name in c(ladder_field_borne_parameters(),
                 ladder_field_through_leaf_area())) {
    got <- ladder_boundary_tangent(patch, name)
    ref <- ladder_boundary_difference(patch, name)

    # Both sides are partial derivatives at one state only if they agree on the
    # value, and a mismatch invalidates the rows below rather than failing them.
    # A few ULP, not zero: the tangent side evaluates the seed's carbon through
    # xad's forward scalar, and on arm64 the compiler fuses that arithmetic
    # differently from the plain double path -- measured 1 ULP (1.1e-16
    # relative) on Apple Silicon where Linux x86-64 reads bit-identical. One
    # state is still one state at that distance; a wrong state is 1e-3 away.
    expect_equal(got$log_density, ref$value[["log_density"]],
                 tolerance = 8 * .Machine$double.eps,
                 label = paste("boundary value,", name))
    expect_equal(got$carbon, ref$value[["carbon"]],
                 tolerance = 8 * .Machine$double.eps,
                 label = paste("birth-size carbon,", name))

    # The seed's height solves live mass equals seed mass, which reads the tissue
    # and allometric parameters and nothing the reduction owns. So k_I's height
    # row is exactly zero for a named reason, and the pair's is not -- which is
    # the same split again, read one quantity earlier.
    reaches_birth_size <- name %in% ladder_field_through_leaf_area()
    if (reaches_birth_size) {
      expect_gt(abs(got$dheight), 0)
    } else {
      expect_identical(got$dheight, 0)
    }

    parts <- c("log_density", "carbon", if (reaches_birth_size) "height")
    for (part in parts) {
      row <- got[[paste0("d", part)]]
      expect_gt(abs(row), 0)
      ladder_report_margin(paste("boundary", part, name),
                           abs(row - ref$row[[part]]) / abs(ref$row[[part]]),
                           1e-3)
    }
  }
})

test_that("the introduction map's whole Jacobian agrees entry by entry", {
  # This fixture's own unit, and the one object it had no reference for. The map is the
  # pre-introduction state and the traits in, the whole widened state out; the
  # forward side seeds one tangent per input column and the reverse side is the
  # transpose the sweep runs. Both go through introduce_over, so the reference
  # traverses the forward function and not the transpose.
  #
  # Formed entirely rather than contracted: a dot product returns one number, hides
  # an error behind a small seed component, and localises to nothing when it fails.
  # At this size the whole object fits.
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29), c(0))
  stand <- ladder_run(p)
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  introduction <- which(diff(widths) != 0)
  times <- vapply(trajectory, function(s) s$time, numeric(1))

  patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), ladder_control())
  for (k in seq_along(introduction)) {
    b <- introduction[[k]]
    pre <- trajectory[[b]]$state
    # The first introduction introduces both species at t = 0; the second is species 1
    # alone, into a patch that already carries nodes -- which is the arrangement a
    # narrowing written as a truncation of the tail fails on.
    who <- if (k == 1L) c(1L, 2L) else 1L
    expect_equal(patch$ode_size, length(pre))

    j <- ladder_introduction_jacobian_tf24(patch, who, pre, times[[b]])
    scale <- max(abs(j$forward))
    expect_gt(scale, 0)
    gap <- max(abs(j$forward - j$reverse)) / scale
    message(sprintf("\n  introduction %d: %d rows x %d columns", k, nrow(j$forward),
                    ncol(j$forward)))
    ladder_report_margin(sprintf("  introduction Jacobian, introduction %d", k),
                         gap, 1e-12)
    for (i in who) patch$introduce_new_node(i, times[[b]])
  }
})

test_that("the state a range resumes from is the one the introduction made", {
  # The state the range above an introduction linearises its first step at is
  # rebuilt rather than recorded -- no record holds it -- and the sweep only
  # length-checks the rebuild. Two things must hold of it and neither was asserted:
  # the rows an introduction does not touch come through unchanged, and the rows it
  # writes are the boundary node the patch was carrying.
  #
  # Tolerance is exactly zero for both: an introduction is a push, not an
  # arithmetic operation on the rows it does not touch.
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29), c(0))
  stand <- ladder_run(p)
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  introduction <- which(diff(widths) != 0)
  times <- vapply(trajectory, function(s) s$time, numeric(1))

  patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), ladder_control())
  for (k in seq_along(introduction)) {
    b <- introduction[[k]]
    pre <- trajectory[[b]]$state
    who <- if (k == 1L) c(1L, 2L) else 1L
    patch$set_ode_state(pre, times[[b]])
    invisible(patch$ode_rates)
    held <- unlist(lapply(who, function(i) patch$species[[i]]$new_node$ode_state))
    before <- vapply(seq_len(length(patch$species)),
                     function(i) patch$species[[i]]$ode_size, numeric(1))
    for (i in who) patch$introduce_new_node(i, times[[b]])
    post <- patch$ode_state

    stride <- length(held) / length(who)
    carried <- integer(0)
    newcomer <- integer(0)
    q <- 0L
    for (i in seq_len(length(patch$species))) {
      carried <- c(carried, q + seq_len(before[[i]]))
      q <- q + before[[i]]
      if (i %in% who) {
        newcomer <- c(newcomer, q + seq_len(stride))
        q <- q + stride
      }
    }
    carried <- c(carried, q + seq_len(length(post) - q))
    expect_identical(post[carried], pre)
    expect_identical(post[newcomer], held)
  }
})

test_that("the initial reserve channel carries a row", {
  # A germination event provisions the seedling: the initial reserve is a trait
  # times the storage capacity, so how well a seed is stocked reads the traits
  # directly. It enters nowhere else, so a dropped channel is an exact zero.
  shared <- ladder_shared("introductions")
  stand <- shared$stand
  result <- shared$gradient
  at <- which(ladder_bare_traits(census_trait_names_tf24(stand)) == "a_st3")
  expect_gt(length(at), 0L)
  peak <- max(abs(result$gradient[, at]))
  message(sprintf("\n  a_st3 row peak: %.6e", peak))
  expect_gt(peak, 0)
})

test_that("the newcomer depends on the state it was introduced into", {
  # The introduction is evaluated at the pre-introduction state: the field is
  # rebuilt and the boundary node placed in it, so the newcomer's own quantities
  # depend on the state before the introduction as well as on the traits. An
  # implementation written as though the newcomer has no predecessor drops that
  # half entirely.
  #
  # Perturbing the soil at the introduction time is what exposes it, because a
  # boundary node taking its potentials from a cache has no visible dependence on
  # moisture at all -- and establishment is the most moisture-sensitive event in
  # the life cycle.
  # The channel is isolated by perturbing the soil in the pre-introduction state
  # itself and introducing from each: comparing two whole runs would move the
  # answer whatever the newcomer reads, because two runs differ everywhere.
  stand <- ladder_shared("introductions")$stand
  # The fixture's own parameters, which the founding patch below is rebuilt from.
  # Building them runs nothing; the stand above is where the run is.
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29, 0.94), c(0, 0.57))

  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  introduction <- which(diff(widths) != 0)
  b <- introduction[[2]]                      # an introduction with a populated field
  pre <- trajectory[[b]]$state
  at <- trajectory[[b]]$time

  # The pre-introduction width is two founding nodes, one per species, which is
  # what the recorded state at this introduction is a state of.
  founding_patch <- function() {
    patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), ladder_control())
    for (i in 1:2) {
      patch$set_time(0)
      patch$compute_environment()
      patch$introduce_new_node(i, 0)
    }
    patch
  }
  newcomer_density <- function(state) {
    patch <- founding_patch()
    patch$set_ode_state(state, at)
    invisible(patch$ode_rates)
    ladder_boundary_density(patch)
  }

  # The soil block begins after every species' nodes, so its offset is read off the
  # patch the recorded state belongs to -- the founding pair, not the stand the run
  # ended at, which has three more nodes and would index past the state entirely.
  template <- founding_patch()
  expect_equal(template$ode_size, length(pre))
  moisture <- sum(vapply(seq_len(length(template$species)),
                         function(i) template$species[[i]]$ode_size,
                         numeric(1))) + 1L

  base <- newcomer_density(pre)
  # The top layer, and it has to be the top one: a plant at birth size is
  # shallow-rooted, so the deepest layer of the column does not reach it and a
  # perturbation there moves the newcomer by exactly zero. That is a property of
  # the seedling rather than of the channel, and a probe placed there would report
  # the channel missing.
  bumped <- pre
  bumped[[moisture]] <- bumped[[moisture]] * 1.05
  moved <- newcomer_density(bumped)

  gap <- max(abs(moved - base) / abs(base))
  message(sprintf("\n  a 5%% perturbation of the top layer moves the boundary node by %.3e",
                  gap))
  expect_gt(gap, 1e-9)
})

test_that("the initial condition carries no trait row on these fixtures", {
  # The last range's state adjoint is the initial-condition term, and the sweep
  # discards the final narrowed vector. That is exact here and it is a scope
  # statement rather than a check passing: every fixture in this ladder starts
  # from bare ground, so the first range is the empty patch and y(0) is the soil
  # alone. No trait reaches it, so there is nothing in the discarded vector.
  #
  # What the discard would lose is a run resumed from a populated exported state,
  # where founding cohorts are in y(0) and their provisioning reads a_st3. This
  # ladder has no such fixture, so the term is unrefereed rather than confirmed,
  # and saying so is the honest form: the earlier version of this test compared
  # two different runs and required the answer to move, which two different runs
  # do whatever the sweep keeps.
  shared <- ladder_shared("introductions")
  stand <- shared$stand
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  # Bare ground: the first recorded width is the narrowest the run ever has, and
  # the very first step is already an introduction, so nothing was integrated before the
  # first introduction. There is no cohort in y(0) for a trait to reach.
  expect_identical(widths[[1]], min(widths))
  expect_gt(widths[[2]], widths[[1]])
  message(sprintf("\n  bare ground: first width %.0f, introduction at step 1",
                  widths[[1]]))
})

test_that("seed height carries a row rather than a declaration", {
  # The seed height solves an implicit condition on the strategy. Preparation runs
  # in double and resolves it before the traits are differentiable inputs, so a
  # strategy that receives the result carries nothing through this channel -- and
  # the loss is exact on a tangent and on a sweep alike, because both inherit the
  # same declaration. That is why it was a declared bias for as long as it was one.
  #
  # It is now derived where the newborn's state is written, so the eight
  # parameters reaching birth size carry a row and the check is a measurement.
  # Every one of them must be a column here; omega is the one whose column is the
  # channel and nothing else, so it is asserted non-zero rather than merely
  # present.
  eight <- ladder_birth_size_parameters()
  expect_length(eight, 8L)

  shared <- ladder_shared("introductions")
  stand <- shared$stand
  columns <- ladder_bare_traits(census_trait_names_tf24(stand))
  expect_length(setdiff(eight, columns), 0L)

  result <- shared$gradient
  at <- which(ladder_bare_traits(colnames(result$gradient)) == "omega")
  expect_length(at, length(stand$patch$species))
  expect_gt(max(abs(result$gradient[, at, drop = FALSE])), 0)
})

test_that("a newcomer's leaf area reaches the census through a field built without it", {
  # The fourth boundary channel, and the only one of the five whose dropped-channel
  # signature is a SHORTFALL rather than an exact zero -- so it is the one that
  # needs a reference rather than a zero test. The newcomer is placed in a field
  # built without it, so its own leaf area reaches the census through a field it
  # does not contribute to, which is a distinct path from every cohort's.
  #
  # The probe is a field-borne parameter, compared column by column against a
  # tangent of the same trajectory. On the introductions fixture every introduction is
  # on the path, so this is also where the introduction boundary's own residual is
  # measured rather than inferred.
  shared <- ladder_shared("introductions")
  stand <- shared$stand
  result <- shared$gradient
  columns <- colnames(result$gradient)

  # The reference has to be a derivative of the same trajectory before its
  # disagreement means anything.
  ladder_report_margin("the replay reaches this stand's own census",
                       ladder_replay_floor(stand, result$value), 1e-9)

  # The census's own reading of the traits, so the trajectory term can be separated
  # from it. For the allometric constants the two very nearly cancel, and the
  # column's relative error is then an amplified view of the trajectory term's.
  own <- do.call(rbind, ladder_census_trait_direct_tf24(stand))
  dimnames(own) <- list(rownames(result$gradient), columns)

  # The field-borne pair, and every parameter reaching the census through the seed's
  # own size. The birth-size class is here because a row missing where a newborn's
  # state is written moved all eight of them together, on this fixture, while the
  # field-borne columns held -- and only two of the eight were being compared, under
  # a bound widened for them. The class is what an introduction puts on the path, so the
  # class is what an introduction's rung asks for.
  wanted <- c("1.k_I", "2.k_I",
              paste0("1.", ladder_birth_size_parameters()),
              paste0("2.", ladder_birth_size_parameters()))
  wanted <- wanted[wanted %in% columns]
  expect_length(wanted, 18L)

  for (name in wanted) {
    reference <- ladder_trajectory_tangent(
      stand, ladder_trait_direction(columns, name))$tangent
    # A shortfall, not a zero: the channel is live in both, and what is asked is
    # whether the sweep's is the same size as the reference's.
    expect_gt(max(abs(result$gradient[, name])), 0)
    r <- ladder_column_residual(result$gradient[, name], reference,
                                result$gradient[, name] - own[, name])
    message(sprintf(
      "    %-8s amplification %6.2f, residual against the trajectory term %.2e, smallest row %.1e of the column",
      name, r$amplification, r$per_row_trajectory, r$share))
    ladder_report_margin(paste("column", name), r$per_row,
                         ladder_trajectory_agreement())
  }
})


test_that("the census's sensitivity to a range's starting state is refereed", {
  # A range's first step runs from a state no record holds -- the widened one,
  # after the introduction and before the step. The three probes here are the only
  # way to reach it: one names the state, one differentiates the whole remaining
  # trajectory from it exactly, and one replays that trajectory in plain double so
  # the first can be differenced.
  #
  # This is the forward counterpart of what the sweep carries and discards: the
  # adjoint it holds after the last introduction is d(census)/d(this state). None
  # of the reverse pass is on this path, so the two share no code.
  stand <- ladder_stand_introductions_short()
  base <- ladder_range_base_state(stand, 0L)

  # Segment 0 is the state no step reached, which on this coordinate holds the
  # environment and no cohort -- so its width is the soil column's.
  expect_gt(length(base), 0L)
  message(sprintf("\n  range 0 starts from %d components", length(base)))

  # An unseeded direction differentiates nothing, which pins that the tangent is
  # carried by the seed rather than by the replay.
  none <- ladder_census_initial_state_tangent_tf24(stand, numeric(length(base)),
                                                   0L)
  expect_true(all(none$tangent == 0))
  expect_true(all(is.finite(none$value)))

  # The replay has to reach the same census the tangent reports, or the two are
  # differentiating different functions.
  #
  # ⚠️ THE BUDGET IS 1e-9 AND IT IS THE SAME CLAIM THIS FILE MAKES ONE TEST ABOVE,
  # where it reads 9.78e-13 against 1e-9. It was 1e-12 here, fitted to a reading of
  # 1.24e-15, and that reading was luck rather than a property: the two routes are
  # one function at two scalars, so with a zero seed they differ only by whatever a
  # nested root-find's iterate count does with a 1e-16 input difference -- which
  # nothing on either route controls, and which reshuffles whenever the curve the
  # find walks changes shape. Measured on that curve's two forms: 1.24e-15 with the
  # conductivity on its own table, 4.4e-13 with it read as the integral's slope,
  # 1.36e-12 with the integral read as a quintic. The tangent-against-difference
  # checks below -- which are what this file is for -- are unmoved by all three.
  replayed <- ladder_census_initial_state_replay_tf24(stand, base, 0L)
  ladder_report_margin("the initial-state replay reaches the run's census",
                       max(abs(replayed - none$value) / abs(none$value)), 1e-9)

  for (i in seq_along(base)) {
    d <- replace(numeric(length(base)), i, 1)
    tangent <- ladder_census_initial_state_tangent_tf24(stand, d, 0L)$tangent
    along <- function(rel) {
      h <- max(abs(base[[i]]), 1) * rel
      (ladder_census_initial_state_replay_tf24(stand, base + h * d, 0L) -
         ladder_census_initial_state_replay_tf24(stand, base - h * d, 0L)) /
        (2 * h)
    }
    coarse <- along(1e-5)
    fine <- along(1e-6)
    scale <- max(abs(fine))
    if (scale == 0) {
      # ⚠️ NOT A SKIP, AND IT WAS ONE. The difference says this component reaches
      # no census metric, and the tangent has to say the same: a live tangent
      # beside a dead difference is a channel one route has and the other does
      # not, which is the disagreement this file exists to find. As a `next` it
      # asserted nothing at all, so a wrong non-zero tangent passed -- and the
      # count is how it hid: the file went 104 passes to 103 when one more
      # component joined this branch, with a message rather than a failure.
      # ⚠️ THE BOUND IS THE CENSUS'S RESOLUTION, NOT ZERO, and the difference
      # between the two readings is what this measures. Components 6 to 9 come
      # back EXACTLY 0.00e+00 -- structurally dead, no channel at all. Component 4
      # comes back 2.2e-11 of the census: a live channel whose response is below
      # what a replayed difference can carry, so the difference reads a hard zero
      # and lands here. Holding it to 1e-12 would report the census's own last
      # bits as a disagreement between the two routes.
      ladder_report_margin(
        sprintf("d(census)/d(range 0 state %d) is zero on both routes", i),
        max(abs(tangent)) / max(abs(none$value)), 1e-9)
      next
    }
    floor <- max(max(abs(coarse - fine)) / scale, 4 * .Machine$double.eps)
    ladder_report_margin(sprintf("d(census)/d(range 0 state %d)", i),
                         max(abs(tangent - fine)) / scale, 10 * floor)
  }
})
