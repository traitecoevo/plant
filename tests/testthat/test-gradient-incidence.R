# The instruments the gradient's scope questions are answered with.
#
# The leaf classifies its operating point by the branch taken and the next plant
# overwrites it, and a clamp that severs a row leaves a number indistinguishable
# from a true zero. Neither is recoverable after a run from anything else, so
# these counters are the only route to "how often" -- and "how often" is what
# decides whether a refused regime is a corner or most of the run.

# One entry per driver, built by ladder_driver_stand() so that a regime this file
# names and a regime the parity file names are the same regime. The ORDER inside
# an entry is the point: a sweep adds to the
# tallies a run leaves, so the forward reading is taken before anything sweeps
# and kept, rather than re-taken from a second run of the same recipe. The
# gradient is taken on first request and kept beside the stand, so a block
# needing both pays for one.
#
# ⚠️ NOTHING MAY CLEAR A SHARED STAND'S DIAGNOSTICS. clear_diagnostics() resets
# every counter, the clamps and the curvature margin together, so a block that
# cleared one to check the clearing left every later reader asserting against
# zero -- which is a check that passes whatever the model does. The readings
# below are copies, so a caller cannot reach the counters at all.
incidence_cache <- new.env(parent = emptyenv())
incidence_key <- function(rain, lifetime, k_I, introductions) {
  paste(rain, lifetime, k_I, introductions)
}

incidence_built <- function(rain, lifetime, k_I = 0.5, introductions = NULL) {
  key <- incidence_key(rain, lifetime, k_I, introductions)
  if (is.null(incidence_cache[[key]])) {
    scm <- ladder_driver_stand(rain, lifetime, k_I, introductions = introductions)
    incidence_cache[[key]] <- list(
      scm = scm,
      kinds = stats::setNames(census_operating_point_counts_tf24(scm)[[1]],
                              census_operating_point_names_tf24()),
      clamps = stats::setNames(census_clamp_counts_tf24(scm)[[1]],
                               census_clamp_names_tf24()))
  }
  incidence_cache[[key]]
}

# The same entry with the sweep run, which is what the differentiated tallies and
# the curvature margin are read off.
incidence_swept <- function(rain, lifetime, k_I = 0.5, introductions = NULL) {
  key <- incidence_key(rain, lifetime, k_I, introductions)
  entry <- incidence_built(rain, lifetime, k_I, introductions)
  if (is.null(entry$gradient)) {
    entry$gradient <- stand_gradient(entry$scm, hyperpar = NULL)
    entry$swept <-
      stats::setNames(census_clamp_counts_differentiated_tf24(entry$scm)[[1]],
                      census_clamp_names_tf24())
    entry$curvature <- census_curvature_margin_tf24(entry$scm)[[1]]
    incidence_cache[[key]] <- entry
  }
  entry
}

incidence_stand <- function(rain, lifetime, k_I = 0.5, introductions = NULL) {
  incidence_built(rain, lifetime, k_I, introductions)$scm
}

# The classification tally as the RUN left it, never as a swept object reports it.
incidence_of <- function(rain, lifetime, k_I = 0.5, introductions = NULL) {
  incidence_built(rain, lifetime, k_I, introductions)$kinds
}

test_that("the classification tally is the route to a regime's incidence", {
  # A wet stand never leaves the branch the gradient answers for, which is what
  # makes it the fixture every other rung uses -- and is why incidence measured
  # on one says nothing about a dry one.
  wet <- incidence_of(2.0, 5)
  expect_gt(wet[["interior"]], 0)
  expect_equal(sum(wet[names(wet) != "interior"]), 0)

  # The tally is cleared and re-accumulated per run rather than carried, or a
  # second measurement would read the first one's states as well. On a stand of
  # its own: clearing is destructive and every other block here reads a shared
  # one. A patch lifetime of 1 is the shortest stand that carries cohorts, so
  # this costs a second.
  scm <- ladder_driver_stand(2.0, 1)
  counts <- census_operating_point_counts_tf24(scm)[[1]]
  expect_gt(sum(counts), 0)
  census_clear_diagnostics_tf24(scm)
  expect_equal(sum(census_operating_point_counts_tf24(scm)[[1]]), 0)
})

test_that("a dry stand is interior throughout, and answers", {
  # The number this exists to produce. It was first taken while this driver's
  # gradient was refused outright, to say how much answering the pinned branch
  # would buy. On phylloptim 0.16 the stem was held at psi_crit, so 54% of this
  # stand's solves were pinned there (`boundary-crit`). psi_crit is no longer a
  # bound and the stem follows the soil past it, so every solve is interior. It
  # is not recoverable from a run afterwards -- the classification is
  # overwritten by the next individual.
  n <- incidence_of(0.25, 10)
  expect_gt(n[["interior"]], 0)
  expect_equal(sum(n[names(n) != "interior"]), 0)

  # And it answers. Counted before it is read, because all(is.finite(x)) is
  # TRUE of an empty vector.
  g <- incidence_swept(0.25, 10)$gradient
  expect_false(any(stand_gradient_refused(g)))
  expect_null(g$refusal[[1]])
  expect_gt(length(g$gradient[[1]]), 0)
  expect_true(all(is.finite(g$gradient[[1]])))
})

test_that("the light floor is counted on both paths, and binds at neither shipped value", {
  # Where this clamp binds, a cohort's radiation stops depending on any other
  # cohort's height. Both halves are asserted, because a counter that never fires
  # and a counter that always fires are equally uninformative.
  #
  # ⚠️ THE FLOOR IS TWO SITES, NOT ONE. The
  # crown site floors each quadrature point of the mean-light integrand; the
  # radiation site floors the mean those points make. The crown weights sum to
  # one, so a mean of floored points cannot fall below the floor: on the shipped
  # shading model only the crown site fires, and reading the radiation site
  # alone would report no severance at all.
  # Scoped to the two light sites: the other clamps have their own regimes and
  # their own block below, and several of them do bind at shipped values.
  nm <- census_clamp_names_tf24()
  light <- match(c("light_floor", "light_floor_crown"), nm)
  shipped <- incidence_built(2.0, 5, k_I = 0.5)$clamps
  expect_true(all(shipped[light] == 0))
  # Non-vacuity, and it is not decorative: this tally is shared and clearing it
  # anywhere would make the line above pass whatever the light field did.
  expect_gt(shipped[["rooting_depth"]], 0)

  # k_I is a free parameter a gradient-driven search walks, and walking it up is
  # what walks the field into the floor. Eighty times the shipped value.
  #
  # Twenty introductions rather than the schedule's eighty-eight: every claim
  # below is about which sites fire and in what order, and the counts reported
  # are shares of this stand's own solves. Measured, the two readings differ in
  # nothing this block asserts and the thinned one costs 18 s against 258.
  entry <- incidence_built(2.0, 5, k_I = 40, introductions = 20L)
  walked <- entry$scm
  fired <- entry$clamps
  solves <- sum(entry$kinds)
  # Counts, not a share of solves: the crown site is counted once per quadrature
  # point of the mean-light integrand, so solves is not its denominator and
  # dividing by it reported 348%.
  message(sprintf("  at k_I = 40, over %.0f solves: %s %.0f, %s %.0f",
                  solves, nm[light[1]], fired[light[1]],
                  nm[light[2]], fired[light[2]]))
  expect_gt(fired[[light[[2]]]], 0)
  # The mean of floored points is at or above the floor, so the radiation site
  # does not bind under mean light.
  expect_equal(fired[[light[[1]]]], 0)

  # The forward model keeps running, and the light floor is not what stops the
  # gradient: below the floor the census is not a function of light at all, so
  # the row is exactly zero for the model as evaluated rather than a row
  # withheld.
  expect_gt(stand_census(walked)[[1]], 0)
  g <- incidence_swept(2.0, 5, k_I = 40, introductions = 20L)$gradient

  # ⚠️ THE FLOOR'S ROW IS A DECLARED ZERO AND NOT A REFUSAL, WHICH IS THE
  # WHOLE POINT OF THIS BLOCK. Below the floor the census is not a function of
  # light at all, so the honest row is the zero -- and this stand answers, so
  # that zero is what the gradient carries rather than what a refusal covered.
  # The sweep left what a double holds on this k_I until the storage pool's
  # charge and drain form came back. Counted before it is read, because
  # all(is.finite(x)) is TRUE of an empty vector.
  expect_false(any(stand_gradient_refused(g)))
  expect_null(g$refusal[[1]])
  expect_gt(length(g$gradient[[1]]), 0)
  expect_true(all(is.finite(g$gradient[[1]])))

  # And the severance is readable rather than silent, which is the whole basis on
  # which the zero is declared instead of refused. The forward tally cannot stand
  # in for this: it counts every solve, where the sweep visits only the recorded
  # steps.
  swept <- incidence_swept(2.0, 5, k_I = 40, introductions = 20L)$swept
  message(sprintf("  the sweep's own severances: %s",
                  paste(sprintf("%s %.0f", nm[light], swept[light]),
                        collapse = "  ")))
  expect_gt(swept[[light[[2]]]], 0)
  expect_equal(swept[[light[[1]]]], 0)
  # A counter a rebind drops reports zero however often the clamp fires, which is
  # what this one did before the storage was shared. And the sweep visits the
  # recorded steps rather than every solve, so its tally is the smaller.
  expect_lt(swept[[light[[2]]]], fired[[light[[2]]]])

  # Named from the enum, so a site cannot be counted under its neighbour's name.
  expect_true(all(c("light_floor", "light_floor_crown") %in% nm))
  expect_length(fired, length(nm))
  expect_length(swept, length(nm))
})

# Every clamp site, classified by what its incidence says rather than by what
# reading the code suggests. The classification lives here rather than in a
# document because a document cannot fail.
#
# ⚠️ THE TEST IS THREE-WAY AND THE FIRST CASE SPLITS IN TWO, which reading the
# code does not reveal:
#
#   the model's own zero   the clamp IS a modelling statement, so the census
#                          genuinely does not depend on what it masks. Declared,
#                          and NOT a candidate for removal -- rooting_depth is
#                          this, and it is the largest severance in the model.
#   a guard's zero         a numerical floor, and the census is bit-identical
#                          either side of it, so the row is exactly zero for the
#                          model as evaluated. Declared, and a candidate for
#                          removal by changing the FORWARD model.
#   never binds            counted, and reported as never having fired, which is
#                          the only thing that separates a guard that held from
#                          one nothing reached.
clamp_class <- list(
  rooting_depth          = "model",
  light_floor            = "guard",
  light_floor_crown      = "guard",
  soil_moisture_floor    = "forward-only",
  soil_potential_ceiling = "forward-only",
  soil_conductivity      = "forward-only",
  soil_positivity        = "never",
  rainfall               = "never",
  infiltration           = "never",
  # `supply_kink_step_off` is deliberately absent: there is no note_clamp site for
  # it. At a coincident collar the supply kernels return not-a-number, and nothing
  # steps four kink tolerances off it -- the marginal profit falls back to a
  # central difference at the same collar, and a row built on a NaN supply
  # derivative is differenced instead of read. Neither is a clamp, so neither is
  # counted here. A site that cannot fire reporting "never" is the one reading
  # this list must not produce.
  # The leaf model's one, the stomatal corner E_s >= 0, which is a different
  # facility: the leaf solves in double on both paths, so it keeps ONE tally and
  # the forward share is the total less the delta measured across
  # record_leaf_outputs. phylloptim measures it active at most solved points on a
  # broad grid, mostly dark, hot or dry leaves. A stand never reaches it, the
  # shaded and clamped drivers' shut leaves included: a shut leaf sits at its wet
  # bound, where E_s is zero by construction rather than by the clamp.
  stomatal_flux          = "never"
)

test_that("every clamp site is classified, and by a measured incidence", {
  nm <- census_clamp_names_tf24()
  # A site with no classification is the drift this list exists to prevent.
  expect_setequal(nm, names(clamp_class))

  wet_s <- incidence_swept(2.0, 5)$swept
  dry_s <- incidence_swept(0.10, 5)$swept
  dry_f <- incidence_built(0.10, 5)$clamps

  message(sprintf("  %-24s %10s %10s", "site", "wet", "drought"))
  for (s in nm) {
    message(sprintf("  %-24s %10.0f %10.0f  [%s]", s, wet_s[[s]], dry_s[[s]],
                    clamp_class[[s]]))
  }

  # Never means never, on both paths and both drivers. A site in this class that
  # starts firing is a regime the suite has never seen — including the leaf's
  # stomatal corner, which a stand does not reach.
  for (s in names(clamp_class)[unlist(clamp_class) == "never"]) {
    expect_equal(wet_s[[s]], 0, label = paste(s, "fired on the wet driver"))
    expect_equal(dry_s[[s]], 0, label = paste(s, "fired on the drought driver"))
    expect_equal(dry_f[[s]], 0, label = paste(s, "fired forward on drought"))
  }

  # Forward-only means the sweep never met it, which is a weaker claim than
  # never and has to be kept separate: the guard is reachable, so a longer run
  # or a finer schedule could put it on a recorded step.
  for (s in names(clamp_class)[unlist(clamp_class) == "forward-only"]) {
    expect_equal(dry_s[[s]], 0,
                 label = paste(s, "reached the differentiated path"))
  }

  # The model's own zero binds on EVERY driver including the control, which is
  # what says it is the model rather than a corner.
  expect_gt(wet_s[["rooting_depth"]], 0)
  expect_gt(dry_s[["rooting_depth"]], 0)

  # And a guard's zero does not bind on the control, which is what makes it a
  # guard rather than the model. Which driver DOES reach the light floor is the
  # light-floor block below, on the k_I the site needs rather than on this one.
  expect_equal(wet_s[["light_floor"]], 0)
  expect_equal(wet_s[["light_floor_crown"]], 0)

  # Both drivers still answer, carrying those declared zeros. A severance that
  # made the gradient wrong would have to show up as a refusal or as a
  # disagreement with a rebuilt difference, and neither is here.
  expect_false(any(stand_gradient_refused(incidence_swept(2.0, 5)$gradient)))
  expect_false(any(stand_gradient_refused(incidence_swept(0.10, 5)$gradient)))
})

test_that("phylloptim's root-vulnerability clamps stay out of reach, with a margin", {
  # These are NOT counted, deliberately, and this is what stands in for a counter.
  #
  # Two clamps in the root curve are thresholds on a layer potential: the curve's
  # argument is clamped into its knot domain at 6.8229 MPa, and the cumulative
  # integral is capped past 7.3132. Both already carry matching derivative kills,
  # so where they bind the row is an honest zero rather than a wrong number -- and
  # the root_b row stays right under the cap by the homogeneity rather than in
  # spite of it. So what is worth asserting is not a count but the DISTANCE, which
  # a counter reading zero cannot report.
  #
  # Instrumenting them would need a phylloptim header edit, hence a reinstall,
  # hence a near-full plant recompile: the cost is the build loop rather than the
  # code. A threshold on a readable quantity does not need it.
  last_knot <- 6.8229
  integral_cap <- 7.3132

  worst <- 0
  for (d in list(list(rain = 2.00, name = "wet"),
                 list(rain = 0.10, name = "drought"),
                 list(rain = 0.05, name = "very-dry"))) {
    scm <- incidence_stand(d$rain, 5)
    e <- scm$patch$environment
    nlayer <- e$get_soil_number_of_depths()
    # ⚠️ THE TRAILING SLOT COUNT IS ASKED FOR, NOT COUNTED. The environment's
    # block is the soil layers followed by its cumulative-flux accumulators, and
    # a window written as `n - 4` skipped four where there are five: at the FIRST
    # recorded state, where the whole vector is the environment's ten entries, it
    # read the first accumulator as the deepest layer's moisture. That slot is
    # exactly 0, `psi_from_soil_moist(0)` returns the 1000 MPa dry cap, and both
    # bounds below then fail on a number no layer ever held -- the worst real
    # potential over these three drivers is 4.34 MPa, on the drought one.
    n_aux <- length(e$get_soil_water_state_cumulative_flux())
    # Over the RECORDED steps, which is the set the sweep visits -- a terminal
    # reading misses a layer that dried and rewetted, and those are the states the
    # clamp would bind in. The forward run keeps them, so this costs no sweep:
    # measured, the three drivers contribute 510, 2503 and 3715 records whether
    # or not a gradient has been taken on them.
    for (r in scm$store_recording()) {
      s <- r$state
      n <- length(s)
      theta <- s[(n - nlayer - n_aux + 1):(n - n_aux)]
      if (length(theta) != nlayer || any(!is.finite(theta))) next
      psi <- vapply(theta, function(x) e$psi_from_soil_moist(x), numeric(1))
      if (any(!is.finite(psi))) next
      worst <- max(worst, max(psi))
    }
  }
  message(sprintf("  worst layer potential over three drivers: %.4f MPa, against %.4f and %.4f",
                  worst, last_knot, integral_cap))
  expect_lt(worst, last_knot)
  expect_lt(worst, integral_cap)

  # Non-vacuity: a run that dried nothing would pass the two bounds above while
  # measuring nothing at all, so the drying has to be real.
  expect_gt(worst, 1.0)
})

test_that("the curvature guard reports how close it came, not only that it held", {
  # A guard that held and a guard nothing reached report the same green, so the
  # distance to the floor is carried out of the run. The floor is a declared
  # Control entry for the same reason: it changes which rows exist, so two
  # gradients taken at different values are gradients of different functions.
  wet <- incidence_swept(2.0, 5)
  margin <- wet$curvature
  floor <- gradient_control(wet$scm)[["gradient_curvature_floor"]]
  message(sprintf("  smallest curvature met: %.4g, against a floor of %.4g (%.0fx)",
                  margin, floor, margin / floor))
  # Non-vacuity both ways: a margin of -1 means the interior branch was never
  # reached, so the reading would say nothing.
  expect_gt(margin, 0)
  expect_gt(margin, floor)
  # And the floor is in the set stand_gradient compares, or two gradients taken
  # at different floors would read as comparable.
  expect_true("gradient_curvature_floor" %in%
                names(gradient_control(wet$scm)))
})
