# Whether the sweep runs at all, and on which channels.
#
# The first block is the ladder's canary: it drives the sweep on the fixture every
# trajectory rung uses, so a sweep that cannot run reports here as well as
# wherever else it is reached.

test_that("the sweep runs on a stand in the declared regime", {
  stand <- ladder_stand_two_by_two()
  ladder_require_regime(stand, "stand")
  g <- stand_gradient(stand)
  # Non-vacuity: a sweep that ran and answered nothing is not a sweep that ran,
  # and every metric refusing is a different event from the sweep failing.
  expect_false(any(stand_gradient_refused(g)))
  expect_gt(sum(is.finite(g$gradient)), 0)
})

test_that("a curvature the collar cannot stand on refuses instead of returning zeros", {
  # The profit row is the envelope theorem's and reads no curvature; every collar
  # response is a quotient over it. So when the curvature is unusable the leaf
  # leaves the water rows off the tape and carries the values on -- and a row left
  # off the tape comes back an exact zero, which is the one thing a caller must
  # not be able to read as an answer.
  #
  # ⚠️ THE GRAIN IS THE WHOLE EVALUATION, NOT THE METRIC, and that is forced
  # rather than chosen: the missing row is an intermediate of one recording that
  # spans six stages and every cohort in them, so no seed carries a component the
  # loss could be attributed to. The refusal therefore travels on the one channel
  # the sweep already has, and a seed that never reads an uptake output pays for it
  # too. What is discharged here is that the loss is REPORTED; separating the two
  # output kinds is not.
  #
  # Injected, because nothing reaches this guard by being driven there -- see
  # ladder_patch_fold for the sweep that establishes that.
  folded <- ladder_patch_fold()
  nm <- ladder_rate_names(folded)
  n <- folded$ode_size

  size_only <- c("height", "mass_heartwood")
  water_coupled <- c("environment_1", "environment_9")

  for (label in c(size_only, water_coupled)) {
    got <- ladder_rhs_adjoint_tf24(folded, replace(numeric(n), match(label, nm), 1))
    expect_true(got$refused,
                label = paste("a seed answered from rows that do not exist:", label))
    expect_true(grepl("curvature", got$reason, fixed = TRUE))
    # Not-a-number rather than the zeros the tape would otherwise hand back, so a
    # caller ignoring the flag still cannot read a severance as an answer. This is
    # the whole reason the flag is not enough on its own.
    expect_true(all(is.na(got$state)))
    expect_true(all(is.na(got$trait)))
  }

  # The control: at the shipped floor the same fixture answers both kinds, so this
  # is measuring the guard rather than the fixture.
  plain <- ladder_patch_one()
  for (label in c(size_only, water_coupled)) {
    got <- ladder_rhs_adjoint_tf24(plain, replace(numeric(n), match(label, nm), 1))
    expect_false(got$refused, label = paste("the unguarded fixture refused:", label))
    expect_true(all(is.finite(got$state)))
    # Non-vacuity: an answer of zeros is what a severed row and a refused one both
    # look like, so the surviving gradient has to be live.
    expect_gt(sum(abs(got$state) > 0), 0)
  }
})

test_that("a fold refuses every metric, and on this census nothing would be spared", {
  # ⚠️ MEASURED, AND IT IS NOT WHAT THE REQUIREMENT ASSUMES. Even a channel that
  # could refuse per metric would save nothing on TF24's census, and the reason is
  # the census rather than the sweep: all three metrics are size moments, but
  # growth reads water, so sweeping backwards gives every metric a non-zero soil
  # adjoint within a step or two of the census. There is no water-independent
  # metric here to be spared.
  #
  # Kept as a check rather than a note because it is a statement about the census
  # that would silently stop being true if a metric were added -- and if one is,
  # this block is where it shows up.
  p <- ladder_parameters(species = "fast", lifetime = 2)
  ctrl <- ladder_control(gradient_curvature_floor = 1e3)
  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.428 * 0.5, 5))
  env$extrinsic_drivers_set_constant("rainfall", 2.0)
  scm <- SCM("TF24", "TF24_Env")(p, env, empty_events(), ctrl)
  scm$run()
  g <- stand_gradient(scm)

  # Every metric refuses, each by the curvature guard rather than by anything
  # else -- so the refusal is the one under test and it reached the metric
  # boundary as data rather than as an exception.
  refused <- stand_gradient_refused(g)
  expect_true(all(refused))
  for (m in rownames(g$gradient)) {
    # The whole row goes, which is what makes one reason enough for it.
    expect_true(all(is.na(g$gradient[m, ])))
    expect_true(grepl("curvature", g$refusal[[m]]$reason, fixed = TRUE))
    # Located, which a refusal carried as data could easily fail to be.
    expect_gte(g$refusal[[m]]$species, 0)
  }
  message(sprintf("  every one of %d metrics refuses at a fold: %s",
                  nrow(g$gradient),
                  paste(rownames(g$gradient), collapse = ", ")))
})

test_that("refusal is metric-level and not per-parameter", {
  # A sum has no defined value with an undefined term, so one refused operating
  # point anywhere in one metric's sweep makes that metric's whole gradient
  # undefined -- not the cohort's column, not the parameter's entry. Metrics are
  # independent of each other, so refusal is metric-level and no wider.
  #
  # A partly-populated gradient vector is the shape that reads as an answer.
  stand <- ladder_stand_two_by_two()
  result <- stand_gradient(stand)
  for (m in rownames(result$gradient)) {
    row <- result$gradient[m, ]
    expect_false(any(is.na(row)) && !all(is.na(row)),
                 label = paste("metric refused for some parameters only:", m))
  }
})

test_that("a shut point answers, and the branch it took is what decides", {
  # Dried to half its moisture the one-cohort patch reaches a genuine hydraulic
  # shutdown, and the shut branches answer for it.
  #
  # The classification is read off the KIND rather than a residual: the marginal
  # profit returns a hard sentinel zero in a no-flow state, and no residual test
  # can tell that from stationarity. It has to come from the exit the solve took --
  # that is what routes this state away from an envelope step it has no stationary
  # point for.
  interior <- ladder_patch_one()
  shutdown <- ladder_patch_shutdown()

  # Both answer, and the shut one is not answered by accident: it is a different
  # derivation, reached because the leaf reports the branch.
  expect_silent(invisible(ladder_block_value_tf24(interior, 1L)))
  expect_silent(invisible(ladder_block_value_tf24(shutdown, 1L)))
})

test_that("a shut point's rows are numbers, and every one of them is finite", {
  # The old form of this check asserted that no vector was produced at all, on
  # the grounds that a refused row coming back as zeros is indistinguishable from
  # an answer. That hazard has not gone away -- it has moved: the rows here ARE
  # numbers now, so what has to hold is that they are finite and that the zeros
  # among them are the structural ones.
  #
  # Where a refusal is still reachable, `test-gradient-parity.R` carries the
  # absence-never-a-number check against a driver that genuinely refuses.
  shutdown <- ladder_patch_shutdown()
  n <- shutdown$ode_size
  seed <- ladder_seeds(n, scale = ladder_block_scale(shutdown$ode_rates))
  got <- ladder_rhs_adjoint_tf24(shutdown, seed)$state
  expect_true(is.numeric(got))
  expect_true(all(is.finite(got)))
  # Non-vacuity: a shut leaf still moves the stand, so a row of zeros would mean
  # the state had been answered for by not being read.
  expect_gt(sum(abs(got) > 0), 0)
})

test_that("both output kinds answer at a shut point, which they did not before", {
  # The two output kinds refuse independently: the profit row survives every
  # degeneracy except a jump of the argmax and an undefined objective, and the
  # uptake row is the one that ceases to exist.
  #
  # At a shut point the question is now moot in the direction that matters --
  # BOTH answer. The rows are a substituted constant plus the bound's movement,
  # and neither output kind is missing. So this measures that all four seeds come
  # back rather than that all four are refused together.
  #
  # ⚠️ The independence requirement is not discharged HERE, and no longer needs
  # to be: the block above discharges it on a fixture where the guard fires.
  # What this one says is narrower and still worth saying -- at a shut point the
  # question is moot, because both kinds answer.
  shutdown <- ladder_patch_shutdown()
  nm <- ladder_rate_names(shutdown)
  n <- shutdown$ode_size
  kinds <- list("size only, height" = "height",
                "size only, heartwood" = "mass_heartwood",
                "water-coupled, a layer" = "environment_1",
                "water-coupled, total uptake" = "environment_9")
  for (label in names(kinds)) {
    seed <- replace(numeric(n), match(kinds[[label]], nm), 1)
    got <- ladder_rhs_adjoint_tf24(shutdown, seed)$state
    message(sprintf("  %-28s answered", label))
    expect_true(all(is.finite(got)))
  }
})
