# The fixtures the reverse-sweep checks run on, their regime assertions, and the
# references every check reports its margin against.
#
# Five things a check has to carry before it counts, and they are why this file
# exists rather than each check holding its own literals:
#
#   path disjointness  every reference written here is in R, over quantities the
#                      forward model publishes. It traverses the forward
#                      reduction, which it must, and no transpose.
#   fault injection    test-gradient-ladder-injection.R corrupts one side of a
#                      comparison the suite already makes and requires the
#                      comparison to reject it, and ladder_report_margin()
#                      records how much room the check had. A fault detected at
#                      three times the tolerance is a check about to stop working.
#   non-vacuity        ladder_expect_moves() zeroes the channel a check exists
#                      for and requires the answer to change.
#   adversariality     the fixtures keep every distinguishing quantity distinct:
#                      no equal heights, no small-integer ratios, two species
#                      differing by whole factors, seeds from a fixed generator.
#   tolerance          ladder_forward_floor() measures a fixture's own arithmetic
#                      floor and every tolerance is a multiple of it.
#
# There are two kinds of fixture and the difference is not cosmetic. A **patch**
# fixture is constructed: its state is written directly, so it can be put in the
# regime the checks need and crossed by hand. A **stand** fixture is run, so its
# state is whatever the trajectory reached and the constructed conditions are not
# available. The block-level checks take a patch; only the checks that need a
# trajectory take a stand.

# Chosen so that no two of the per-species reduction parameters are equal or in a
# small-integer ratio between the species. A reduction sum that collapses the
# species into one scalar, or a scatter that reaches the wrong species, then
# moves a number.
ladder_traits <- function() {
  list(
    fast = c(lma = 0.0825, hmat = 5.13,  k_I = 0.5,  a_l1 = 5.44, a_l2 = 0.306),
    slow = c(lma = 0.1074, hmat = 11.71, k_I = 0.31, a_l1 = 3.61, a_l2 = 0.401))
}

ladder_parameters <- function(species = "fast", birth_rate = NULL,
                              lifetime = 2, k_I = NULL) {
  traits <- ladder_traits()
  if (is.null(birth_rate)) {
    birth_rate <- c(fast = 1.10, slow = 0.83)[species]
  }
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  for (i in seq_along(species)) {
    tr <- traits[[species[[i]]]]
    # A shading regime is one extinction coefficient for the whole stand, not one
    # per species: the two carry different values by default, so leaving one alone
    # would make the regime a difference between the species as well as a change
    # in the light.
    if (!is.null(k_I)) {
      tr[["k_I"]] <- k_I
    }
    p <- add_strategies(p, trait_matrix(unname(tr), names(tr)),
                        hyperpar = TF24_hyperpar,
                        birth_rate = list(unname(birth_rate[[i]])))
  }
  p
}

ladder_control <- function(...) {
  Control(node_density_in_birth_date = TRUE, ...)
}

# The patch fixtures, whose state is written rather than reached.

# One species, one cohort: the smallest stand that still carries endogenous
# feedback, because it shades itself and draws on its own soil layers. Both
# reduction transposes, the retention factor and the field's slope channel are
# live here. Two cohorts buy accumulation, not feedback, and keeping them apart
# is what makes a failure at four nodes localisable.
ladder_patch_one <- function(parameters = NULL) {
  ladder_patch(species = "fast",
               heights = list(4.73),
               log_densities = list(-0.39),
               parameters = parameters)
}

# The one-cohort patch dried until its leaf stops having an interior optimum.
#
# No injected fault: at half this fixture's moisture the operating point is a
# genuine hydraulic shutdown, which is one of the five operating-point kinds and
# one the boundary refuses by name. That makes a refusal reachable from a
# constructed patch in about a second, where a refusing RUN is not available at
# all -- the trajectory rungs' stands stay interior by construction, which is
# what their regime asserts.
#
# `factor` is how far down; the default is the first value at which the refusal
# fires, so a check on it is measuring the branch and not a deep extrapolation.
# The one-cohort patch with the curvature floor raised above every curvature the
# model produces, so the guard that separates the two output kinds fires at a
# state the model actually reaches.
#
# It is an INJECTION and it has to be, which is the finding rather than a
# shortcut. Sweeping 5625 solved leaf states found 1351 interior ones, every one
# strictly concave and the smallest magnitude 0.0623: the folds are real -- at
# stem_c 0.4, 40% of the feasible interval is non-concave -- but the solve pins
# rather than landing on one, and at a pin the curvature is never computed. So no
# fixture reaches this guard by being driven there, and raising the floor is the
# only way to establish that it does what its name says.
ladder_patch_fold <- function(floor = 1e3) {
  ladder_patch(species = "fast",
               heights = list(4.73),
               log_densities = list(-0.39),
               control = ladder_control(gradient_curvature_floor = floor))
}

ladder_patch_shutdown <- function(factor = 0.5) {
  heights <- list(4.73)
  densities <- list(-0.39)
  ladder_condition(
    ladder_patch(species = "fast", heights = heights, log_densities = densities),
    heights, densities,
    moisture = c(0.192, 0.271, 0.233, 0.317, 0.208) * factor)
}

# Three cohorts of one species, with the birth dates chosen so that the first and
# third grid points carry EQUAL trapezium weights.
#
# That equality is what makes a permutation available. The reductions are
# quadratures over the birth-date abscissa, so exchanging the states of two grid
# points changes the sum unless their weights agree -- and with weights
# (x2-x1)/2, (x3-x1)/2, (x4-x2)/2, (x4-x3)/2 over the ascending grid the first and
# third agree exactly when x2-x1 = x4-x2. The boundary node is the fourth point,
# born at the clock, so x4 is the fixture's own time and x2 is its midpoint with
# x1.
#
# `swap` exchanges the two states. Every shared object then sees the same two
# plants in the opposite order, with the same field and the same soil, which is the
# one thing a re-run cannot vary.
# The one-cohort patch given a root system that puts the collar on the soil.
#
# The uniform-drying direction is a near-symmetry only when the collar follows a
# uniform shift of the soil potentials -- measured, 0.973 of it here against
# 0.698 at the shipped traits -- and when root mass is spread over the layers, so
# the top layer's own column does not already carry the tracking. Both are
# properties of the root network and neither is reached by a trajectory: every
# run stand measures 1.1x, and a run with these traits ends at 1.0x.
#
# `conductance` scales root mass and divides root respiration and turnover by the
# same factor, so the root CARBON budget is bit-for-bit the shipped one and only
# the hydraulics differ. Net production is 2.09 against the shipped 1.25.
# The default is the first value at which the amplification passes fifteen; the
# leaf refuses between 100 and 130.
ladder_patch_uniform_drying <- function(conductance = 60, spread = 0.99) {
  p <- ladder_parameters("fast")
  strategies <- p$strategies
  for (i in seq_along(strategies)) {
    pars <- strategies[[i]]$pars
    pars$a_r1 <- pars$a_r1 * conductance
    pars$r_r  <- pars$r_r / conductance
    pars$k_r  <- pars$k_r / conductance
    pars$root_depth_shape_eta <- spread
    strategies[[i]]$pars <- pars
  }
  p$strategies <- strategies
  heights <- list(4.73)
  densities <- list(-0.39)
  ladder_condition(
    ladder_patch(species = "fast", heights = heights,
                 log_densities = densities, parameters = p),
    heights, densities)
}

ladder_patch_permutable <- function(swap = FALSE, time = 1.37) {
  births <- c(0.17, (0.17 + time) / 2, 1.03)
  stopifnot(births[[2]] < births[[3]], births[[3]] < time)
  heights <- c(4.73, 6.29, 3.11)
  densities <- c(-0.39, -1.67, -2.71)
  order <- if (swap) c(3L, 2L, 1L) else c(1L, 2L, 3L)
  ladder_patch(species = "fast",
               heights = list(heights[order]),
               log_densities = list(densities[order]),
               birth_dates = list(births),
               time = time)
}

# Each node's own rates, keyed by the height it carries rather than by where it
# sits. A permutation moves a state to a different position, so a reader indexing
# by position could not tell the two runs apart.
ladder_node_rate_table <- function(x) {
  patch <- ladder_as_patch(x)
  invisible(patch$ode_rates)
  out <- list()
  for (i in seq_len(length(patch$species))) {
    sp <- patch$species[[i]]
    for (j in seq_along(sp$nodes)) {
      node <- sp$nodes[[j]]
      out[[length(out) + 1L]] <- list(species = i, position = j,
                                      height = node$height,
                                      rates = node$individual$ode_rates)
    }
  }
  out
}

# How far apart a permutation leaves the shared field, which is this fixture's own
# floor rather than a tolerance anyone chose.
#
# Exchanging two grid points of equal weight leaves the reduction's SUM invariant
# and its summation ORDER reversed, so the knot values move in their last bits and
# bit-identity is not available at this level. Measured on the permutable fixture
# the values agree to 6e-10 and the knot heights exactly; the slopes carry a larger
# relative figure only where the slope is near zero, so the value channel is what a
# floor is taken from.
#
# Anything a carried quantity would do is orders above this.
ladder_permutation_floor <- function(a, b) {
  ka <- ladder_field_knots(a)
  kb <- ladder_field_knots(b)
  stopifnot(identical(ka[, 1], kb[, 1]))
  scale <- max(abs(ka[, 2]), .Machine$double.xmin)
  max(max(abs(ka[, 2] - kb[, 2])) / scale, 4 * .Machine$double.eps)
}

# The field the recorded step reads: one row per knot, height then value then
# slope. Off the published spline, which is where the run leaves it -- the
# reduction that builds the field is an intermediate of the stage recording and
# has no transpose of its own, so what a knot serves here is the forward check
# that permuting the nodes leaves it alone.
ladder_field_knots <- function(x) {
  ladder_as_patch(x)$environment$light_availability$state
}

# The state a range's first step ran from, which is what a first-range reference
# indexes its direction against.
#
# The run keeps it and nothing needs to rebuild it: store_trajectory()'s first
# record is the first range's base, and a later range's base is the introduction
# record that widened into it. Checked bit for bit against the C++ accessor this
# replaced over 45 ranges on four stands.
ladder_range_base_state <- function(x, range = 0L) {
  records <- x$store_trajectory()
  if (range == 0L) {
    return(records[[1L]]$state)
  }
  Filter(function(r) isTRUE(r$introduction), records)[[range]]$state
}

# Two species, two cohorts each. Two species catch a reduction sum that collapses
# them into one scalar; two cohorts of one species catch a trait accumulated per
# cohort rather than per species. Four nodes catch both.
#
# `cross` puts one pair out of birth-date order in height, which is the state the
# forward model runs and a guard written on height order refuses. It is on by
# default because folding the condition into every run is the point of it;
# passing FALSE is how a failure is localised to the crossing.
ladder_patch_two_by_two <- function(cross = TRUE, parameters = NULL) {
  heights <- if (cross) list(c(4.73, 6.29), c(5.41, 2.57))
             else       list(c(6.29, 4.73), c(5.41, 2.57))
  ladder_patch(species = c("fast", "slow"),
               heights = heights,
               log_densities = list(c(-0.39, -1.67), c(-1.33, -2.71)),
               parameters = parameters)
}

# A live patch with its state written to the values handed in. Nodes are
# introduced oldest-first, which is the order the state vector holds them in, and
# each at its own time: a node takes its birth date from the clock, so
# introducing them all at one instant leaves the abscissa tied and every question
# about birth-date order unanswerable.
ladder_patch <- function(species, heights, log_densities,
                         relative_reserve = 0.12, time = 1.37,
                         birth_dates = NULL, parameters = NULL,
                         control = NULL) {
  p <- if (is.null(parameters)) ladder_parameters(species) else parameters
  ctrl <- if (is.null(control)) ladder_control() else control
  patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), ctrl)
  if (is.null(birth_dates)) {
    offset <- c(0, 0.17, 0.11)
    birth_dates <- lapply(seq_along(heights), function(i) {
      offset[[i]] + seq(0, by = 0.29, length.out = length(heights[[i]]))
    })
  }
  events <- do.call(rbind, lapply(seq_along(species), function(i) {
    data.frame(species = i, time = birth_dates[[i]])
  }))
  for (k in order(events$time)) {
    # The birth date is the inserter's own argument. Routed through the clock
    # instead, every node of a species shares a date and the abscissa the
    # birth-date coordinate integrates over is tied.
    patch$set_time(events$time[[k]])
    patch$compute_environment()
    patch$introduce_new_node(events$species[[k]], events$time[[k]])
  }
  ladder_condition(patch, heights, log_densities, relative_reserve, time)
}

# Put the stand where the checks need it.
#
# The reserve placement is the vacuity guard and is the least obvious of these.
# The gate is centred near the bottom of the reserve range with a width of a
# tenth, so at comfortable reserves its slope is orders below its slope in the
# band. A stand left where a run put it damps every channel running through
# growth, and a check then passes because the signal is small rather than because
# the code is right.
ladder_condition <- function(patch, heights, log_densities,
                             relative_reserve = 0.12, time = 1.37,
                             moisture = c(0.192, 0.271, 0.233, 0.317, 0.208)) {
  state <- patch$ode_state
  at <- 0L
  for (i in seq_len(length(patch$species))) {
    sp <- patch$species[[i]]
    stride <- sp$nodes[[1]]$ode_size
    names_i <- sp$nodes[[1]]$ode_names
    i_height <- match("height", names_i)
    i_storage <- match("storage", names_i)
    i_density <- match("log_density", names_i)
    stopifnot(length(heights[[i]]) == length(sp$nodes),
              length(log_densities[[i]]) == length(sp$nodes))
    for (j in seq_along(sp$nodes)) {
      base <- at + (j - 1L) * stride
      h <- heights[[i]][[j]]
      state[[base + i_height]] <- h
      state[[base + i_density]] <- log_densities[[i]][[j]]
      state[[base + i_storage]] <-
        relative_reserve * ladder_storage_capacity(patch, i, h)
    }
    at <- at + sp$ode_size
  }
  # The layers are given distinct moistures for the same reason the heights are
  # distinct: a uniform column makes a scatter that reaches the wrong layer, and
  # a retention factor applied at the wrong one, invisible.
  n_layer <- patch$environment$get_soil_number_of_depths()
  stopifnot(length(moisture) >= n_layer)
  state[at + seq_len(n_layer)] <- moisture[seq_len(n_layer)]
  patch$set_ode_state(state, time)
  # The auxiliaries the regime is read off are written by a rate evaluation, not
  # by loading a state, so a report taken without this reads the previous one's.
  invisible(patch$ode_rates)
  patch
}

# The stand fixtures, which need a trajectory.

# Two species, two cohorts each, run rather than constructed. The trajectory
# rungs need the steps the adaptive pass resolved, so this cannot be conditioned;
# what it can be is chosen, and the birth rates are chosen to keep every cohort
# out of the absorbing reserve region.
ladder_stand_two_by_two <- function() {
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.63), c(0, 0.41))
  ladder_run(p)
}

# Three introductions in the order species 1, species 2, species 1, so the node
# stride is exercised in both directions. Introducing into species 1 while
# species 2 exists is the arrangement a narrowing implemented as a truncation of
# the tail fails on.
ladder_stand_introductions <- function() {
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29, 0.94), c(0, 0.57))
  ladder_run(p)
}

# The introductions fixture with every species' establishment constant scaled, so
# a recruit can be put where the establishment probability is stiff.
#
# The probability is P^2/(P^2 + k^2) times a decay in patch age, with
# k = a_d0 * a_0, so its derivative in P peaks at P = k/sqrt(3) -- ratio 1 below.
# At the shipped constant the fixture's recruits sit at ratio 36 and 21, which is
# the flat region where the check would prove nothing; the scale is what moves
# them onto the peak.
# The same introductions over a quarter of the run, for the checks whose assertion is
# an exact identity rather than a measured margin.
#
# Composition over steps is associative or it is not; a sweep is repeatable or it is
# not; a permutation of the metric order changes the numbers or it does not. None of
# those is a statement about run length, and all three are asserted with
# expect_identical at tolerance zero -- so a shorter fixture tests the same claim and
# re-blesses nothing. What it must keep is the SHAPE: five introductions in the same
# order, so the node stride is still exercised in both directions and the split still
# has interior steps either side of an introduction to cut at.
#
# The schedule is compressed rather than truncated for that reason, and the dates
# stay mutually non-commensurate.
ladder_stand_introductions_short <- function() {
  p <- ladder_parameters(c("fast", "slow"), lifetime = 0.45)
  p$node_schedule_times <- list(c(0, 0.13, 0.31), c(0, 0.22))
  ladder_run(p)
}

# Many ranges, and a short run. The sweep opens a range at every introduction, so
# the range count is what the range machinery and the insertion transpose are
# exercised at -- and in the DEFAULT schedule that count is derived from the patch
# lifetime, which is why the two read as one thing. They are separable: written
# introduction times put 61 ranges in a run of less than half a year.
#
# Measured on this fixture: 62 introductions, 62 ranges, 68 accepted steps, and a
# sweep in 1.7 s. The four-node stand every other trajectory rung uses has six.
#
# ⚠️ WHAT THIS DOES NOT REACH IS THE STEP COUNT. Sixty-eight against the ~3,400 a
# century fixture accepts, so nothing here says the recording's footprint or the
# step loop behaves at production length. That gap is recorded in NEWS.md under
# Known issues and is not closed by this fixture.
#
# It is outside the declared regime on one assertion and DELIBERATELY SO: at this
# cohort count the heights pile up and stop being non-commensurate. That
# separation is what makes a swapped scatter target visible to a localisation
# reference, and a bit-identity consults no reference at all -- so what the
# fixture gives up is exactly what the check it carries does not use. Do not call
# ladder_require_regime() on it, and do not use it for anything that reads a
# margin.
#
# The dates are a golden-ratio sequence rather than an even spread: consecutive
# terms of one are never in a small-integer ratio, which is the separation
# `ladder_non_commensurate` asks for, and the sequence is deterministic so the
# fixture is the same on every machine.
ladder_stand_many_ranges <- function(n = 61L, lifetime = 0.45) {
  phi <- (sqrt(5) - 1) / 2
  times <- sort(round(((seq_len(n) * phi) %% 1) * lifetime, 5))
  p <- ladder_parameters(c("fast", "slow"), lifetime = lifetime)
  # The second species takes every other date, so the two never widen together and
  # the node stride is exercised in both directions at every width.
  p$node_schedule_times <- list(c(0, times), c(0, times[c(TRUE, FALSE)]))
  ladder_run(p)
}

# A run that starts from a patch already carrying cohorts, so the recording holds
# seventeen steps before the state first widens.
#
# Every other stand introduces at the initial time, which puts the first introduction
# at the first recorded step and leaves the lowest range with no step in it. A
# walk that ran that range and a walk that started above it are then the same
# walk. Here they are not: half the run sits below the first introduction, and the
# census depends on the state there through three cohorts rather than through the
# soil alone.
#
# The birth dates are supplied rather than taken from make_initial_state, which
# gives every seeded node the same one; the birth-date coordinate integrates over
# that abscissa and refuses a tie. They are distinct, negative and mutually
# non-commensurate for the reason every other fixture's are.
ladder_stand_resumed <- function() {
  p <- ladder_parameters(c("fast", "slow"), lifetime = 0.45)
  p$node_schedule_times <- list(c(0, 0.13, 0.31), c(0, 0.22))
  state <- make_initial_state(p,
                              heights = list(c(1.71, 0.93), c(1.24)),
                              log_densities = list(c(-0.39, -1.67), c(-1.03)),
                              env = Environment("TF24"), ctrl = ladder_control())
  state$node_times <- list(c(-0.37, -0.11), c(-0.23))
  ladder_run(set_initial_state(p, state))
}

# Which rows of a recording are introductions, and how many of its rows are steps.
#
# A recording holds one row per instruction, and an introduction is an instruction: it
# shares its time with the row below it and no step reached it. So the number of
# steps is not the number of rows, and the row an introduction happens at is a row of
# its own rather than the gap between two.
ladder_introduction_rows <- function(trajectory) {
  which(vapply(trajectory, function(s) s$introduction, logical(1)))
}

# The first row is where the run started and no step reached it either, so it comes
# off with the introductions.
ladder_step_count <- function(trajectory) {
  length(trajectory) - length(ladder_introduction_rows(trajectory)) - 1L
}

ladder_stand_marginal_recruit <- function(scale = 30, lifetime = 0.45) {
  p <- ladder_parameters(c("fast", "slow"), lifetime = lifetime)
  strategies <- p$strategies
  for (i in seq_along(strategies)) {
    pars <- strategies[[i]]$pars
    pars$a_d0 <- pars$a_d0 * scale
    strategies[[i]]$pars <- pars
  }
  p$strategies <- strategies
  p$node_schedule_times <- list(c(0, 0.13, 0.31), c(0, 0.22))
  ladder_run(p)
}

# Where each species' recruit sits on the establishment probability's stiff band,
# as a multiple of the potential at which its derivative peaks. One is the peak.
ladder_recruit_stiffness <- function(x) {
  patch <- ladder_as_patch(x)
  vapply(seq_len(length(patch$species)), function(i) {
    node <- patch$species[[i]]$new_node
    production <- node$individual$aux("net_mass_production_dt")
    k <- ladder_strategy_parameter(patch, i, "a_d0") *
      ladder_area_leaf(patch, i, node$height)
    production / (k / sqrt(3))
  }, numeric(1))
}

ladder_recruit_establishment <- function(x) {
  patch <- ladder_as_patch(x)
  vapply(seq_len(length(patch$species)), function(i) {
    patch$species[[i]]$new_node$individual$establishment_probability(
      patch$environment)
  }, numeric(1))
}

# The smallest stand that shows the allometric-constant disagreement, and its
# control: one species, four tenths of a year, with and without a second cohort.
#
# It exists so work on that defect has a loop measured in seconds. The pair runs in
# about six seconds each against forty-three for the two-species fixture the
# measurements were first taken on, and it separates the defect by a factor of six
# thousand -- 1.4e-03 against 2.2e-07 on the trajectory term. A second cohort
# cannot be had without an introduction after t = 0, so the control differs by that
# introduction as well as by the cohort; every other check here has excluded the
# introduction, which is what makes the pair a usable control rather than a
# confound.
ladder_stand_allometric_probe <- function(two_cohorts = TRUE) {
  p <- ladder_parameters("fast", lifetime = 0.4)
  p$node_schedule_times <- list(if (two_cohorts) c(0, 0.29) else 0)
  ladder_run(p)
}

# One species introduced once, at t = 0, so the recorded trajectory never widens
# and every state it holds can be loaded back into one patch. That is what lets a
# per-stage check run at the states a run actually visited rather than only at a
# state written by hand.
ladder_stand_trajectory <- function(lifetime = 2) {
  p <- ladder_parameters("fast", lifetime = lifetime)
  p$node_schedule_times <- list(0)
  ladder_run(p)
}

ladder_stand_one <- function() {
  p <- ladder_parameters("fast")
  p$node_schedule_times <- list(c(0, 0.63))
  ladder_run(p)
}

ladder_run <- function(p, ctrl = NULL, env = NULL) {
  if (is.null(ctrl)) ctrl <- ladder_control()
  if (is.null(env)) env <- Environment("TF24")
  run_scm(p, env, ctrl, collect = FALSE)
}

# The water a run is driven with. `rain = NULL` is the default environment, so a
# fixture that does not ask for a regime is unchanged by this existing.
#
# The knots are per year rather than per run, so the seasonality a driver realises
# does not depend on the lifetime it is run for -- at a fixed knot count a longer
# run would realise a slower cycle and the regime would stop being the regime.
ladder_environment <- function(rain = NULL, amplitude = 0, lifetime = 2) {
  env <- Environment("TF24")
  if (is.null(rain)) {
    return(env)
  }
  env$set_soil_water_state(rep(0.428 * 0.5, 5))
  if (amplitude == 0) {
    env$extrinsic_drivers_set_constant("rainfall", rain)
  } else {
    t <- seq(0, lifetime,
             length.out = max(200L, as.integer(ceiling(lifetime * 48))))
    env$extrinsic_drivers_set_variable(
      "rainfall", t, pmax(rain * (1 + amplitude * sin(2 * pi * t)), 0))
  }
  env
}

# A stand driven by one of the regimes below, on the birth-date coordinate the
# reverse sweep transposes. One builder, because the incidence file and the
# parity file ask different questions of the SAME five drivers, and two copies of
# a recipe are two regimes as soon as one of them is edited.
#
# `introductions` thins the default schedule to that many, spread over it. The
# regime is the water and the light; the schedule is the caller's, and what it
# decides is COST rather than regime: the default schedule derives its
# introduction count from the lifetime, so a stand that wants the stiff regime
# pays for cohorts it is not asking about. Measured on the clamped driver at a
# lifetime of 5: twenty introductions record 3267 steps and run and sweep in
# 17.8 s where its eighty-eight record far more and take 258, with both light
# sites firing on both paths and the same range refusal either way.
ladder_driver_stand <- function(rain, lifetime, k_I = 0.5, amplitude = 0,
                                introductions = NULL) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  tr <- c(lma = 0.0825, hmat = 5.13, k_I = k_I, a_l1 = 5.44, a_l2 = 0.306)
  p <- add_strategies(p, trait_matrix(unname(tr), names(tr)),
                      hyperpar = TF24_hyperpar, birth_rate = list(1.10))
  if (!is.null(introductions)) {
    full <- p$node_schedule_times[[1]]
    p$node_schedule_times <-
      list(full[round(seq(1, length(full), length.out = introductions))])
  }
  ctrl <- Control()
  ctrl$node_density_in_birth_date <- TRUE
  scm <- SCM("TF24", "TF24_Env")(p, ladder_environment(rain, amplitude, lifetime),
                                 empty_events(), ctrl)
  # The tallies accumulate, so a stand starts its own count rather than the
  # count of whatever ran before it in this process.
  census_clear_diagnostics_tf24(scm)
  scm$run()
  scm
}

# The regimes a reference is captured over, and the reason each is here rather
# than a sweep over rainfall: each one is chosen to reach a different kind of
# operating point, and which kinds a regime actually reached is recorded beside
# the numbers rather than assumed. The first three move the water and are what
# reach an interior, a pinned and a shut point; the last two move the light and
# are what reach the floor the clamp counter sees.
ladder_reference_regimes <- function() {
  list(
    list(name = "wet",      rain = 2.00, amplitude = 0),
    list(name = "drought",  rain = 0.10, amplitude = 0),
    list(name = "seasonal", rain = 1.00, amplitude = 1.0),
    list(name = "shaded",   rain = 2.00, amplitude = 0, k_I = 20),
    list(name = "clamped",  rain = 2.00, amplitude = 0, k_I = 40)
  )
}

# Reading a stand.

# A stand hands out a snapshot of its patch; a patch is itself. Reading through
# one accessor is what lets the references below serve both fixtures.
ladder_as_patch <- function(x) {
  if (is.null(x$species)) x$patch else x
}

# Every node of every species, boundary node included. The boundary node is the
# reduction's lower grid point and is not ODE state, so a reader that stops at
# the introduced nodes is reading a different distribution from the one the model
# integrates.
ladder_nodes <- function(x) {
  patch <- ladder_as_patch(x)
  out <- list()
  for (i in seq_len(length(patch$species))) {
    sp <- patch$species[[i]]
    nodes <- c(list(sp$new_node), rev(sp$nodes))
    out[[i]] <- data.frame(
      species = i,
      boundary = c(TRUE, rep(FALSE, length(nodes) - 1L)),
      birth_date = vapply(nodes, function(n) n$introduction_time, numeric(1)),
      height = vapply(nodes, function(n) n$height, numeric(1)),
      log_density = vapply(nodes, function(n) n$log_density, numeric(1)))
  }
  do.call(rbind, out)
}

# Each species' inflow boundary density, as the patch currently holds it. The
# boundary node is not ODE state, so which evaluation of the condition a patch is
# carrying is a property of what was last done to it, and that is what the two
# readers below are for.
ladder_boundary_density <- function(x) {
  patch <- ladder_as_patch(x)
  vapply(seq_len(length(patch$species)),
         function(i) patch$species[[i]]$new_node$log_density, numeric(1))
}

# The inflow condition at each of its two evaluations, from one recorded state.
#
# A stage evaluates it twice. `in_field` is the one a state load produces, taken
# in the field with every species' boundary interval left off -- the field is then
# rebuilt including it, so this is the value the reductions were built on.
# `in_uptake` is the one a rate evaluation produces, taken in that rebuilt field;
# it is the value the water aggregation reads, the value an introduced node
# inherits, and the value a census reports.
#
# They are the same function at different arguments, so nothing about either number
# says which one a caller is holding.
ladder_boundary_evaluations <- function(patch, state, time) {
  patch$set_ode_state(state, time)
  in_field <- ladder_boundary_density(patch)
  invisible(patch$ode_rates)
  list(in_field = in_field, in_uptake = ladder_boundary_density(patch))
}

# A gradient's columns with their species index stripped. A check written about a
# strategy's parameter finds it in either species' column through this; a check
# about one species' column names it in full.
ladder_bare_traits <- function(x) sub("^[0-9]+\\.", "", x)

ladder_strategy_parameter <- function(x, species, name) {
  ladder_as_patch(x)$species[[species]]$nodes[[1]]$individual$strategy$pars[[name]]
}

# Leaf area from height, the allometry written out here rather than read back
# from the model: a reference that asked the model for its own integrand would
# share the code it is refereeing.
ladder_area_leaf <- function(x, species, height) {
  a_l1 <- ladder_strategy_parameter(x, species, "a_l1")
  a_l2 <- ladder_strategy_parameter(x, species, "a_l2")
  (height / a_l1)^(1 / a_l2)
}

# The capacity the reserve gate is built on.
ladder_storage_capacity <- function(x, species, height) {
  par <- function(n) ladder_strategy_parameter(x, species, n)
  eta <- par("eta")
  eta_c <- 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta)
  area_sapwood <- ladder_area_leaf(x, species, height) * par("theta")
  par("a_st1") * area_sapwood * height * eta_c * par("rho")
}

# The trapezium rule, signed, so a grid running the wrong way subtracts exactly
# as the forward reduction's does.
ladder_trapezium <- function(x, y) {
  n <- length(x)
  stopifnot(n >= 2L, length(y) == n)
  sum(diff(x) * (y[-1] + y[-n])) * 0.5
}

# The leaf-area census over a named coordinate, summed over species.
#
# A census is a quadrature of a density, so its weights are gaps in the
# coordinate the state is carried on. This forms the same sum over either choice
# of abscissa, so the two can be compared against the model's own and the answer
# says which axis the model integrated.
ladder_census_leaf_area <- function(x, coordinate = c("birth_date", "height")) {
  coordinate <- match.arg(coordinate)
  nodes <- ladder_nodes(x)
  total <- 0
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i, , drop = FALSE]
    weighted <- exp(part$log_density) * ladder_area_leaf(x, i, part$height)
    axis <- part[[coordinate]]
    ord <- order(axis)
    total <- total + ladder_trapezium(axis[ord], weighted[ord])
  }
  total
}

# The census's direct term for the leaf-area metric: the metric's own reading of
# the traits, at fixed state. It is not a sensitivity of the state at all, so no
# sweep produces it, and it is easy to omit precisely because it is a one-line
# calculation at the final state.
#
# Leaf area is a power of height over the two allometric constants, so the two
# rows are elementary and are written here rather than asked of the model.
#
# One entry per species per parameter, named as the gradient's columns are. A
# species' own row reads only its own cohorts, so summing the species and
# comparing the total against a column would compare a sum of two things against
# one of them.
ladder_census_direct_term <- function(x) {
  nodes <- ladder_nodes(x)
  out <- numeric(0)
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i, , drop = FALSE]
    a_l1 <- ladder_strategy_parameter(x, i, "a_l1")
    a_l2 <- ladder_strategy_parameter(x, i, "a_l2")
    area <- ladder_area_leaf(x, i, part$height)
    d_a_l1 <- -area / (a_l2 * a_l1)
    d_a_l2 <- -area * log(part$height / a_l1) / (a_l2^2)
    density <- exp(part$log_density)
    axis <- part$birth_date
    ord <- order(axis)
    term <- c(a_l1 = ladder_trapezium(axis[ord], (density * d_a_l1)[ord]),
              a_l2 = ladder_trapezium(axis[ord], (density * d_a_l2)[ord]))
    out <- c(out, stats::setNames(term, paste0(i, ".", names(term))))
  }
  out
}

# The rows of one rate evaluation, named. Nodes come species-major in the order
# the state vector holds them and the environment closes it, so a Jacobian row
# index can be read as the rate it belongs to.
ladder_rate_names <- function(x) {
  patch <- ladder_as_patch(x)
  out <- character(0)
  for (i in seq_len(length(patch$species))) {
    for (node in patch$species[[i]]$nodes) {
      out <- c(out, node$ode_names)
    }
  }
  n <- length(patch$ode_state)
  stopifnot(length(out) <= n)
  c(out, paste0("environment_", seq_len(n - length(out))))
}

ladder_aux <- function(x, name) {
  patch <- ladder_as_patch(x)
  out <- numeric(0)
  for (i in seq_len(length(patch$species))) {
    for (node in patch$species[[i]]$nodes) {
      out <- c(out, node$individual$aux(name))
    }
  }
  out
}

ladder_net_production <- function(x) ladder_aux(x, "net_mass_production_dt")

ladder_relative_reserve <- function(x) {
  patch <- ladder_as_patch(x)
  out <- numeric(0)
  for (i in seq_len(length(patch$species))) {
    for (node in patch$species[[i]]$nodes) {
      ind <- node$individual
      storage <- ind$ode_state[[match("storage", ind$ode_names)]]
      capacity <- ladder_storage_capacity(patch, i, node$height)
      out <- c(out, min(max(storage, 0) / capacity, 1))
    }
  }
  out
}

# The gate's slope as a fraction of its slope at the band's centre, so the floor
# is a share of the strongest response rather than an absolute number.
ladder_reserve_gate_slope <- function(x) {
  patch <- ladder_as_patch(x)
  r <- ladder_relative_reserve(patch)
  out <- numeric(0)
  at <- 1L
  for (i in seq_len(length(patch$species))) {
    a_st2 <- ladder_strategy_parameter(patch, i, "a_st2")
    for (node in patch$species[[i]]$nodes) {
      g <- 1 / (1 + exp(-(r[[at]] - a_st2) / 0.1))
      out <- c(out, (g * (1 - g)) / 0.25)
      at <- at + 1L
    }
  }
  out
}

ladder_light_at_reads <- function(x) {
  patch <- ladder_as_patch(x)
  out <- numeric(0)
  for (i in seq_len(length(patch$species))) {
    for (node in patch$species[[i]]$nodes) {
      out <- c(out, patch$environment$get_environment_at_height(node$height))
    }
  }
  out
}

# The third reference: a plain-double difference of the rates.

# The whole right-hand side differenced in each ODE state entry, in plain double.
#
# This is a reference the forward tangent cannot be, and the reason is structural
# rather than a matter of accuracy. The environment holds its integrated state as
# a double store and takes it passively, so a tangent seeded on a soil moisture
# carries nothing: every soil column and every soil rate row of its Jacobian is
# exactly zero. A difference re-runs the forward balance, so it carries them.
#
# It is a legitimate referee here, and only here, because the soil balance has no
# supplied row in it. Where a row IS supplied -- the leaf's -- a difference of the
# step that consumes it returns zero whether the row is right, wrong or absent,
# which is why the block's rows are refereed by the tangent instead.
# `relative` chooses which step, and the choice is a statement about the fixture
# rather than a preference.
#
# The default sizes the step max(|state|, 1) * rel, which is relative for a
# component of order one and absolute below it. On a CONSTRUCTED patch every
# component is either of order one or exactly zero, so that is the derivative at
# the value in question and the figures every blessed tolerance here was measured
# against.
#
# Along a TRAJECTORY it is not. The reserve and both heartwood accumulators run
# from 1e-11 down to 1e-25 at early states, so an absolute step of 1e-6 moves them
# by many orders and the column it returns is a difference over a range the model
# is not locally linear on. Those columns read as enormous Jacobian errors and the
# instrument is what produced them. `relative = TRUE` is the form to use there, at
# the price of a smaller step -- and therefore more round-off -- on every component
# below order one.
#
# `refereeable` marks the columns that carry a relative step under either setting.
# A component that is exactly zero has none, and a caller forming this matrix along
# a trajectory reports which columns those were rather than quoting a residual over
# all of them.
ladder_rhs_state_difference <- function(x, rel = 1e-6, relative = FALSE) {
  patch <- ladder_as_patch(x)
  n <- patch$ode_size
  state <- patch$ode_state
  time <- patch$ode_time
  at <- function(v) { patch$set_ode_state(v, time); patch$ode_rates }
  out <- matrix(0, n, n)
  refereeable <- rep(TRUE, n)
  for (j in seq_len(n)) {
    h <- if (relative) abs(state[[j]]) * rel else max(abs(state[[j]]), 1) * rel
    if (!(abs(state[[j]]) > 0) || state[[j]] + h == state[[j]]) {
      refereeable[[j]] <- FALSE
      if (!(h > 0)) h <- rel
    }
    up <- state; up[[j]] <- state[[j]] + h
    dn <- state; dn[[j]] <- state[[j]] - h
    out[, j] <- (at(up) - at(dn)) / (2 * h)
  }
  # Leave the fixture where it was found, and rebuild the auxiliaries a regime
  # report is read off.
  invisible(at(state))
  attr(out, "refereeable") <- refereeable
  out
}

# The finer difference's own error, measured rather than chosen.
#
# `coarse` is taken at ten times `fine`'s step. A central difference is
# second-order while truncation dominates, so the coarse step's error is about a
# hundred times the fine one's and the gap between them is therefore the COARSE
# error. Dividing by that factor is what turns the gap into an estimate of the
# error of the step actually being used -- taking the gap itself as the tolerance
# leaves it loose by two orders, which would not catch a soil row wrong by ten.
ladder_difference_floor <- function(coarse, fine, order = 2) {
  scale <- max(max(abs(coarse)), max(abs(fine)), .Machine$double.xmin)
  gap <- max(abs(coarse - fine)) / scale
  max(gap / 10^order, 4 * .Machine$double.eps)
}

# The same parameter set with one species' strategy parameter moved.
#
# A round trip through the list is what writes it through: the accessors hand out
# copies, so the nested value has to be reassigned outward at each level.
ladder_perturbed_parameters <- function(species, index, name, value) {
  p <- ladder_parameters(species)
  strategies <- p$strategies
  pars <- strategies[[index]]$pars
  pars[[name]] <- value
  strategies[[index]]$pars <- pars
  p$strategies <- strategies
  p
}

# The whole right-hand side differenced in one strategy parameter by REBUILDING
# the strategy from its parameters, at the state the base fixture holds.
#
# This is the complement of the prepared-strategy difference, and the two are
# unbiased on disjoint sets. Perturbing a prepared strategy cannot reach a trait
# the leaf holds its own copy of, so that difference reads exactly zero there. A
# rebuild runs preparation, so the leaf receives the trait -- at the price of also
# moving the seed height, which the differentiated path imposes to zero.
#
# The price is not paid on the columns this is for. The seed height solves
# mass_live(h) = omega, which reads the allometric and tissue parameters and no
# leaf trait at all. So for exactly the columns the prepared perturbation is blind
# to, a rebuild carries no birth-size channel either, and each scheme is unbiased
# precisely where the other is not.
ladder_rate_difference_rebuilt <- function(base, index, name, rel = 1e-6) {
  state <- base$ode_state
  time <- base$ode_time
  value <- ladder_strategy_parameter(base, index, name)
  h <- max(abs(value), 1) * rel
  at <- function(v) {
    patch <- ladder_patch_two_by_two(
      cross = FALSE,
      parameters = ladder_perturbed_parameters(c("fast", "slow"), index, name, v))
    # The base fixture's own state, imposed rather than reconditioned, so nothing
    # but the parameter differs between the two evaluations.
    patch$set_ode_state(state, time)
    patch$ode_rates
  }
  (at(value + h) - at(value - h)) / (2 * h)
}

# The inflow boundary's own quantities and their row in one registered parameter,
# with the parameter named rather than counted.
ladder_boundary_tangent <- function(patch, name, species = 1L) {
  at <- match(paste0(species, ".", name), ladder_trait_names_tf24(patch))
  stopifnot(!is.na(at))
  ladder_boundary_density_tangent_tf24(patch, at)
}

# The same quantities from a strategy REBUILT at a perturbed parameter, at the
# state the patch holds, differenced. The patch is the one-cohort fixture, which
# is what the rebuild reconstructs.
#
# Rebuilding runs preparation, so the seed's size moves as it does on the
# differentiated path, and the state is imposed rather than reconditioned, so both
# sides are partial derivatives at one state. The free check that they are is that
# they agree on the value, and every caller asserts it before reading a row.
# What a central difference of a bisected root-find can be worth.
#
# height_seed() runs util::uniroot at offspring_production_tol, which is passed as
# BOTH the absolute and the relative tolerance, and returns the midpoint of the
# final bracket. So the bracket is at most tol * (1 + root) wide, each evaluation
# is at most half of that from the root, and a central difference of two of them
# carries the whole width over its own step.
#
# This error RISES as the step falls, so the gap between two steps reads the
# opposite of it and a floor taken that way is blind to the whole of it.
ladder_bisection_bound <- function(patch, species, root, step, row) {
  tol <- patch$parameters$strategies[[species]]$control$offspring_production_tol
  tol * (1 + root) / (2 * step * abs(row))
}

ladder_boundary_difference <- function(patch, name, species = 1L, rel = 1e-4) {
  state <- patch$ode_state
  time <- patch$ode_time
  value <- ladder_strategy_parameter(patch, species, name)
  h <- abs(value) * rel
  at <- function(v) {
    built <- ladder_patch_one(
      ladder_perturbed_parameters("fast", species, name, v))
    built$set_ode_state(state, time)
    built$compute_environment()
    invisible(built$ode_rates)
    node <- built$species[[species]]$new_node
    c(log_density = node$log_density, height = node$height,
      carbon = node$individual$aux("net_mass_production_dt"),
      area_leaf = node$individual$aux("competition_effect"))
  }
  up <- at(value + h)
  dn <- at(value - h)
  list(value = at(value), row = (up - dn) / (2 * h))
}

# Which registered parameters reach the light field only through a cohort's cached
# leaf area, and which reach it as parameters of the reduction itself. The two
# groups fail differently when an auxiliary is derived at the wrong scalar, and
# telling them apart is what localises such a failure to the auxiliary.
ladder_field_through_leaf_area <- function() c("a_l1", "a_l2")
ladder_field_borne_parameters <- function() c("k_I")

# The leaf's own traits, which no difference of the rates can referee.
#
# The leaf holds its own copy of each of these, taken when the strategy was
# prepared, and a rate evaluation does not push them back in. So perturbing the
# strategy's parameter leaves the physiology at the value it was prepared with and
# the difference comes back EXACTLY zero -- not small, zero -- whether the sweep's
# row is right, wrong or absent. The same shape as a supplied row one level down.
#
# Their referee is the forward model rebuilt from its parameters, in
# `test-gradient-ladder-factorisation.R`, which loops over exactly this set. A
# rebuild runs preparation, so the leaf is CONSTRUCTED with the moved trait
# rather than handed it afterwards, and it reaches these columns because the
# reference used everywhere else is silent on them.
ladder_leaf_own_traits <- function() {
  # ⚠️ Maintained, not derived, and nothing structural keeps it current: a trait
  # whose row becomes supplied joins the set this asserts and has to be added
  # here. That has now happened three times -- dark respiration, then the two
  # photosynthetic capacities -- and each time the failure was this assertion
  # rather than anything about the gradient, which is the good case. The general
  # form: a list read off a signature goes silently incomplete when the
  # signature grows.
  # Renamed with the (P50, c) reparameterisation: `c` and `b` are the STEM pair
  # and are `stem_c` and `stem_P50`; `root_b` is `root_P50`; `beta2` carries its
  # curve's name. `g1_TF24` is gone -- the leaf takes no Medlyn slope on this
  # path. Exactly the drift this comment predicts, arriving a fourth time.
  c("stem_c", "stem_P50", "TF24_beta2", "TF24_cost_scale", "a",
    "curv_fact_elec_trans", "curv_fact_colim", "vcmax_25", "jmax_25",
    "R_d_25", "root_c", "root_P50")
}

# The regime a fixture has to sit in.

# Each row is one assertion, whether it holds, and its measured value. The
# fixtures are written rather than reached, so a violated assertion is the model
# having moved under a fixture that did not: ladder_require_regime() fails on it.
#
# `level` is which set applies. A constructed patch can be put anywhere, so it
# carries the whole set. A run trajectory cannot be crossed by hand and cannot
# have its reserves placed, so those two are reported rather than enforced on one
# -- and that gap is itself worth knowing, because it is the part of the regime
# the trajectory rungs do not exercise.
ladder_regime_report <- function(x, level = c("patch", "stand")) {
  level <- match.arg(level)
  patch <- ladder_as_patch(x)
  env <- patch$environment
  nodes <- ladder_nodes(patch)

  n_layer <- env$get_soil_number_of_depths()
  moisture <- env$get_soil_water_state()[seq_len(n_layer)]
  sat <- rep_len(env$soil_moist_sat, n_layer)
  residual <- 1e-2

  production <- ladder_net_production(patch)
  light <- ladder_light_at_reads(patch)

  entry <- function(name, ok, value) {
    data.frame(assertion = name, ok = isTRUE(ok),
               value = paste(signif(value, 6), collapse = " "),
               stringsAsFactors = FALSE)
  }

  report <- rbind(
    entry("soil moisture strictly interior on every layer",
          all(moisture > 1.5 * residual) && all(moisture < 0.99 * sat),
          moisture),
    entry("net production strictly positive at every cohort",
          length(production) > 0 && all(production > 0), min(production)),
    entry("light at every read at least 100 times the floor",
          length(light) > 0 && all(light > 100 * 1e-4), min(light)),
    entry("no two heights equal or in a small-integer ratio",
          ladder_non_commensurate(nodes$height), length(nodes$height)),
    entry("no two layers at the same moisture",
          level != "patch" || ladder_non_commensurate(moisture), moisture),
    entry("both species separated in every reduction parameter",
          ladder_species_separated(patch), ladder_species_ratio(patch)))

  # The vacuity guard, and it is measured at both levels because dropping it at one
  # is how a fixture comes to pass on a damped signal. A constructed patch is put
  # in the band; a run stand lands where the trajectory took it, which is a
  # relative reserve near a half and a gate slope an order below the floor -- so
  # every channel running through growth is damped there by roughly the ratio of
  # the two slopes. Reported at stand level rather than enforced, because moving
  # the stands into the band retunes every margin the suite has recorded.
  reserve <- ladder_relative_reserve(patch)
  gate <- ladder_reserve_gate_slope(patch)
  report <- rbind(
    report,
    entry("relative reserve inside the gate's transition band",
          length(reserve) > 0 && all(reserve > 0.02) && all(reserve < 0.30),
          range(reserve)),
    entry("reserve gate slope above its floor",
          length(gate) > 0 && all(gate > 0.4), min(gate)))
  report$enforced <- TRUE
  if (level == "stand") {
    at <- report$assertion %in%
      c("relative reserve inside the gate's transition band",
        "reserve gate slope above its floor")
    report$enforced[at] <- FALSE
  }
  # Crossing is an adversariality condition rather than a regime one, and the
  # fixtures come in both forms, so the check that needs it asserts it itself.
  report
}

ladder_require_regime <- function(x, level = c("patch", "stand")) {
  report <- ladder_regime_report(x, level)
  unmet <- report[!report$ok & !report$enforced, , drop = FALSE]
  if (nrow(unmet) > 0L) {
    message("  regime measured but not enforced here: ",
            paste(sprintf("%s (%s)", unmet$assertion, unmet$value),
                  collapse = "; "))
  }
  bad <- report[!report$ok & report$enforced, , drop = FALSE]
  # FAILS rather than skips. The fixtures are written, not reached, so leaving a
  # regime is the model moving under a fixture that did not -- which is the one
  # event this suite exists to report. A skip here reads as green.
  testthat::expect_equal(
    nrow(bad), 0L,
    label = paste0("fixture outside its declared regime: ",
                   paste(bad$assertion, collapse = "; ")))
  invisible(report)
}

ladder_non_commensurate <- function(x, tol = 1e-3) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2L) return(TRUE)
  for (i in seq_along(x)) {
    for (j in seq_along(x)) {
      if (i >= j) next
      if (abs(x[[i]] - x[[j]]) < tol * max(abs(x[[i]]), abs(x[[j]]))) return(FALSE)
      ratio <- max(x[[i]], x[[j]]) / min(x[[i]], x[[j]])
      if (any(abs(ratio - 2:4) < 1e-2)) return(FALSE)
    }
  }
  TRUE
}

# Pairs whose birth-date order and height order disagree. Reserve-gated growth
# lets a younger cohort overtake an older one, so a transpose guarding on height
# order refuses exactly the stands the forward model runs correctly.
ladder_crossing_count <- function(nodes) {
  total <- 0L
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i, , drop = FALSE]
    part <- part[order(part$birth_date), , drop = FALSE]
    if (nrow(part) < 2L) next
    total <- total + sum(diff(part$height) > 0)
  }
  total
}

ladder_species_ratio <- function(x) {
  patch <- ladder_as_patch(x)
  if (length(patch$species) < 2L) return(NA_real_)
  vapply(c("k_I", "a_l1", "a_l2"), function(n) {
    ladder_strategy_parameter(patch, 1, n) / ladder_strategy_parameter(patch, 2, n)
  }, numeric(1))
}

ladder_species_separated <- function(x) {
  ratio <- ladder_species_ratio(x)
  if (all(is.na(ratio))) return(TRUE)
  all(abs(log(ratio)) > 0.2)
}

# Block-normalised seeds from a fixed generator. Soil moisture is of order a
# third and heartwood mass is in kilograms, so an unnormalised inner product is a
# test of the largest block alone; and a seed of all ones makes a summation-order
# defect invisible. Signs alternate so a sum that should cancel does.
ladder_seeds <- function(n, scale = NULL, seed = 20240811L) {
  had <- exists(".Random.seed", envir = globalenv())
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  v <- stats::runif(n, 0.3, 1.7) * rep_len(c(1, -1), n)
  if (!is.null(scale)) {
    stopifnot(length(scale) == n)
    v <- v / ladder_block_scale(scale)
  }
  v
}

# The magnitude each entry is normalised by: its own size where it has one, and
# the block's median otherwise, so a zero entry does not divide.
ladder_block_scale <- function(x) {
  nonzero <- abs(x[x != 0 & is.finite(x)])
  m <- if (length(nonzero) > 0) stats::median(nonzero) else 1
  if (!is.finite(m) || m == 0) m <- 1
  ifelse(abs(x) > 0 & is.finite(x), abs(x), m)
}

# Tolerance and margin.

# The phrase the SWEEP's own refusal carries, as distinct from every refusal the
# leaf raises.
#
# A descent is a product of step Jacobians and has no error control, so it can
# pass far outside the range of the answer it returns and come back: on the
# drought stand it reaches 4.99978e+281 and contracts to 3.69e+46 one range
# later, for a gradient of order 1e+03. A driver that answers does so with
# margin, and one that overflows is refused rather than handed the overflow's
# numbers.
#
# ⚠️ ONE STRING, because three files key on it -- the parity gate's known-gap
# list, the incidence file's two regimes, and the reference's per-regime skip.
# Three copies would let a change to the message pass two of them.
ladder_range_refusal <- "left the representable range"

# The fixture's own arithmetic floor: how far apart two evaluations of the same
# forward quantity land when the shared solver has been driven elsewhere in
# between. Every tolerance is a multiple of this rather than a literal, so a
# threshold that starts failing can be argued about rather than only raised.
#
# This is also the permutation check's cheapest form: a quantity that does not
# come back to itself is carrying something between solves.
ladder_forward_floor <- function(x) {
  patch <- ladder_as_patch(x)
  state <- patch$ode_state
  time <- patch$ode_time
  first <- patch$ode_rates
  patch$set_ode_state(state * 1.05, time)
  invisible(patch$ode_rates)
  patch$set_ode_state(state, time)
  second <- patch$ode_rates
  scale <- pmax(abs(first), abs(second))
  drift <- ifelse(scale > 0, abs(first - second) / scale, 0)
  max(max(drift), 4 * .Machine$double.eps)
}

ladder_margin <- function(observed, tolerance) {
  if (tolerance <= 0) return(Inf)
  observed / tolerance
}

# Compare and report how much room the check had. A fault caught at three times
# the tolerance is a check about to stop working, so the number is recorded and
# not only compared.
# The statistic the entry-by-entry Jacobian checks compare, in one place, so a
# corruption injected into one side is measured by the same number the real check
# uses rather than by a second one written to resemble it.
#
# Each cell is scaled by the larger of its two rows. A cell whose true value is
# orders below its row carries no information about the row, and a per-cell measure
# would make the tolerance a statement about that cell's round-off; taking the
# larger of the two rows keeps a row that is near zero on one side only from
# reporting a ratio of round-offs.
ladder_matrix_residual <- function(observed, reference) {
  stopifnot(identical(dim(observed), dim(reference)))
  row_scale <- pmax(apply(abs(reference), 1, max), apply(abs(observed), 1, max),
                    .Machine$double.xmin)
  max(abs(observed - reference) / row_scale)
}

# The other direction of ladder_report_margin: a corruption the check is supposed
# to catch has to be caught, and by how much is the evidence the check works.
#
# A fault detected at three times the tolerance is a check about to stop working,
# so the factor is reported and not only the verdict.
ladder_report_detection <- function(label, observed, tolerance) {
  message(sprintf("  %-52s %9.2e / %9.2e  detected at %8.1fx",
                  label, observed, tolerance,
                  observed / max(tolerance, .Machine$double.xmin)))
  testthat::expect_gt(observed, tolerance)
}

ladder_report_margin <- function(label, observed, tolerance) {
  message(sprintf("  %-52s %9.2e / %9.2e  uses %6.3f of budget",
                  label, observed, tolerance,
                  ladder_margin(observed, tolerance)))
  testthat::expect_lt(observed, tolerance)
}

# Zero the channel a check exists for and require the answer to move. A check on
# the slope channel passes with the whole slope channel missing if the fixture's
# slope adjoints happen to be zero.
ladder_expect_moves <- function(with_channel, without_channel, label,
                                floor = 1e-8) {
  denom <- max(max(abs(with_channel)), .Machine$double.eps)
  moved <- max(abs(with_channel - without_channel)) / denom
  testthat::expect(
    moved > floor,
    sprintf(paste("%s: zeroing the channel moved the answer by only %.3e, so",
                  "the check would pass with that channel absent"),
            label, moved))
  invisible(moved)
}

# One sweep per fixture, for the checks that only read the gradient back.
#
# A sweep of the four-node stand is 43.6 s and building the stand is 0.13 s, so the
# sweep is the whole cost of a trajectory check; twenty of them pay one purely to
# read a column out of it. The sweep is deterministic -- the floor asserts two
# consecutive sweeps of one recording are bit-identical -- so a shared one carries
# the same numbers as a private one.
#
# What a check may NOT take from here is a sweep when the sweep IS its subject.
# Repeatability, the metric permutation, the interior split and
# record-once-sweep-many each have to drive stand_gradient themselves, and they do.
#
# Keyed by fixture name rather than by the object, because two calls to a
# constructor return two objects carrying the same trajectory. Under testthat's
# parallel runner every file is its own process, so the cache never spans files and
# a stand cannot arrive at a check carrying another file's history.
#
# A sweep that raises is a failure of the file that asked for it. DO NOT catch it
# here and skip: a refusal the fixtures cannot reach is not a condition to carry
# machinery for, and a sweep that started raising would turn the trajectory tier
# quiet rather than red.
ladder_shared_cache <- new.env(parent = emptyenv())

ladder_shared <- function(key) {
  if (!exists(key, envir = ladder_shared_cache, inherits = FALSE)) {
    stand <- switch(
      key,
      two_by_two = ladder_stand_two_by_two(),
      introductions = ladder_stand_introductions(),
      marginal_recruit = ladder_stand_marginal_recruit(),
      allometric_probe = ladder_stand_allometric_probe(TRUE),
      stop("ladder_shared: no such fixture: ", key, call. = FALSE))
    assign(key, list(stand = stand, gradient = stand_gradient(stand)),
           envir = ladder_shared_cache)
  }
  get(key, envir = ladder_shared_cache, inherits = FALSE)
}

# Fault injection.

# The declared zero lists.

# Every trait column resolves to exactly one of these, and the list is part of
# the check. An exact zero is the signature of a missing accumulator and never of
# true insensitivity, which makes it indistinguishable from an ecological
# finding, so zeros are classified rather than tolerated.

# A parameter belongs on one of the lists below only when it has no live route at
# all; one whose birth-size channel alone is imposed to zero has a live column and
# belongs on the birth-size list further down.
#
# A declared zero is asserted and not merely tolerated, and each list is a claim
# about *why* its columns are zero that ladder_zero_cause() measures at the rate
# level: the cause names which rows of the right-hand side the parameter moves, so
# a column that is zero for a reason other than the one declared fails as loudly
# as an undeclared zero does.
#
# The lists are per level, because the levels have different zeros. A column zero
# in the cohort block can be live in the census and the reverse, so excusing one
# level's zero at the other is how an absent row gets excused.

# The dry bound of the leaf's feasible interval, and nothing inside it. The bound
# is the lesser of the stem's and the root's critical potentials, and at an
# interior operating point the solve never reaches it: moving one leaves profit,
# the collar and every per-layer draw bit-identical, so the row is exactly zero
# by complementary slackness rather than by omission.
#
# This is a property of the regime and not of the model. At a pinned operating
# point these two carry the whole row -- and where the root's potential wins the
# minimum, that row is exactly minus the unit vector in its own direction. The
# fixtures are pinned interior, so nothing here exercises it.
ladder_zero_at_an_interior_optimum <- function() {
  # ⚠️ EMPTY, AND NOT BECAUSE THE ARGUMENT STOPPED HOLDING. psi_crit and
  # root_psi_crit were here while they were free parameters PASSED to the leaf.
  # They are derived from (P50, c) inside it now and the strategy hands over
  # neither, so there is no column to classify -- the complementary-slackness
  # reasoning above is still true of the bound and no longer describes a trait.
  #
  # Kept as a named empty list rather than deleted: the reasoning is what a pinned
  # operating point would need, and a curve whose critical potential becomes free
  # again belongs here.
  character(0)
}

# The accessory cost of a seed. It reaches two rates and no others -- offspring
# production, and the survival-weighted offspring the census does not read -- so
# the column is exactly zero on this metric set and would be live on a fitness
# functional.
#
# omega was here while its second route, through birth size, was imposed to zero.
# The seed height now solves its own condition, so that route carries a row: omega
# sets the height a newborn is given, and through it the storage the height scales
# and every resource row the newborn's own physiology writes.
ladder_zero_outside_the_metric_support <- function() {
  "a_f3"
}

# The rates a parameter outside the metric support is allowed to move. Naming
# them is what makes the class a measurement: a third rate would mean the
# census's silence about the column is wrong.
ladder_reproductive_rates <- function() {
  c("fecundity", "offspring_produced_survival_weighted")
}

# Reach the census through the introduction boundary and through nothing else.
# The block is one individual's physiology and reads none of them: the
# establishment probability, its decay in patch age and the initial reserve are
# evaluated once per introduction, at the seed size, in the field the newcomer is
# placed in. So an exact zero in the BLOCK is correct here and an exact zero in
# the census is not -- which is the whole content of the class, and is why the
# two levels carry different declared lists.
#
# The other half of the claim is asserted where it can be: the floor's census
# classification puts all three in the live class, and the introductions file's single-channel
# probes take a_st3 and recruitment_decay by name.
ladder_zero_outside_the_cohort_block <- function() {
  c("a_d0", "a_st3", "recruitment_decay")
}

ladder_zero_by_construction <- function() {
  c(ladder_zero_at_an_interior_optimum(),
    ladder_zero_outside_the_metric_support())
}

# Which cause a declared zero is claimed to have. Named so a check can measure
# the claim rather than only record it.
ladder_zero_cause <- function(name) {
  bare <- ladder_bare_traits(name)
  if (bare %in% ladder_zero_at_an_interior_optimum()) {
    "a limit the operating point is away from: moves no rate at all"
  } else if (bare %in% ladder_zero_outside_the_metric_support()) {
    "outside the metric support: reaches no rate the census reads"
  } else {
    NA_character_
  }
}

# The trajectory reference: a tangent run of the same solve, stepped at the sizes
# the run recorded. `direction` carries one weight per gradient column, so a unit
# vector asks for one column and a mixed one for a contraction.
#
# It returns the metrics the replay itself reached as well as the derivative,
# because a reference whose value disagrees with the model is a reference to a
# different trajectory and every derivative taken from it is that trajectory's.
ladder_trajectory_tangent <- function(stand, direction) {
  ladder_trajectory_tangent_tf24(stand, direction)
}

# One column's disagreement with a reference, per metric row and normalised per
# row.
#
# Normalising a whole column by its largest row is what a single number invites and
# it hides the small rows: leaf area and above-ground mass are of order one while
# stem area is of order 1e-4, so a stem-area row wrong by three per cent reads as
# 1e-6 of the column's peak. Residuals go per block for exactly this reason, and a
# column is a block here.
#
# `trajectory` is the column with the census's direct term removed. Where the total
# is a near-cancellation of the two terms -- which it is for the allometric
# constants -- the relative error against the total is an amplified view of the
# error in the trajectory term, and the second number is the unamplified one.

# The floor under the per-row scale, as a fraction of the column's own peak.
#
# Per-row normalisation is what keeps the small rows visible, and dividing a cell
# by ITSELF all the way down does the opposite. The two failure modes need
# separating rather than trading, which is what a floor does.
#
# Measured on the introductions fixture over the eighteen shortlisted columns: every
# stem-area cell sits between 7e-07 and 1.1 of its column's peak except one, the
# second species' root constant, whose stem-area entry is 8.4e-09 against a peak of
# 3.2e-01 -- 2.6e-08 of it. The sweep and the tangent carry about 1e-10 of absolute
# disagreement on that row whatever the column, so dividing by 8.4e-09 reports 9e-03
# where the same disagreement reports 5e-06 two rows up.
#
# 1e-5 is the largest fraction at which every other cell in that table is still
# refereed against its own value, and it keeps the case the per-row rule exists for:
# a three per cent error on a stem-area cell of ordinary size is 3e-02 of its own
# scale, a hundred times the bound it is held to.
ladder_column_floor <- function() 1e-5

ladder_column_residual <- function(got, reference, trajectory = NULL) {
  floor <- ladder_column_floor()
  scale <- pmax(abs(reference), max(abs(reference)) * floor)
  out <- list(per_row = max(abs(got - reference) / scale),
              # Which row set the figure, and how far below its column that row
              # sits. A cell the floor had to rescue is one this column does not
              # carry, and saying so is what stops the number reading as a defect.
              at = which.max(abs(got - reference) / scale),
              share = min(abs(reference) / max(max(abs(reference)),
                                               .Machine$double.xmin)))
  if (!is.null(trajectory)) {
    out$per_row_trajectory <-
      max(abs(got - reference) / pmax(abs(trajectory),
                                      max(abs(trajectory)) * floor))
    out$amplification <- max(abs(trajectory) / scale)
  }
  out
}

# The reference a tangent cannot be.

# A differenced column moves one registered parameter and no other.
#
# A trait set through the hyperparameter function moves several: lma carries leaf
# turnover and leaf dark respiration with it, and each of those is a separate
# gradient column. A difference taken that way is the derivative along a
# direction in parameter space, and the sweep's column is the derivative along a
# basis vector. Both are finite, plausible and the same sign, and on lma at four
# years they differ by a factor of three -- so the refusal belongs where the
# perturbation is made, because nothing downstream of it can tell the two apart.
ladder_assert_one_parameter <- function(before, after, name) {
  same <- before == after[names(before)]
  moved <- names(before)[is.na(same) | !same]
  if (identical(moved, name)) {
    return(invisible(moved))
  }
  stop(sprintf("differencing %s moved %d registered parameters (%s)", name,
               length(moved), paste(moved, collapse = ", ")), call. = FALSE)
}

# A whole-run central difference of the census: set one registered parameter, run
# the model, and difference. It shares no code with the sweep or the tangent and,
# which is the point, it inherits none of their declarations. Setting a parameter
# re-runs preparation, so the seed height moves, where both differentiated paths
# hold it still -- and that is the one channel neither of them can referee,
# because a tangent imposes the same equation the sweep does.
#
# It is for the columns a declared zero reaches, not for a sweep: two model runs
# per column against one seeded run for a tangent column.
#
# The parameter is set on the built strategy rather than on the trait vector, and
# the seed-height slope is identical either way -- so this keeps the whole of what
# the reference is for and drops the only thing it cannot referee.
ladder_run_difference <- function(traits, name, rel = 1e-5, lifetime = 0.4,
                                  times = c(0, 0.29), birth_rate = 1.10) {
  build <- function() {
    p <- scm_base_parameters("TF24")
    p$max_patch_lifetime <- lifetime
    p <- add_strategies(p, trait_matrix(unname(traits), names(traits)),
                        hyperpar = TF24_hyperpar,
                        birth_rate = list(birth_rate))
    p$node_schedule_times <- list(times)
    p
  }
  census_at <- function(value) {
    p <- build()
    before <- unlist(p$strategies[[1]]$pars)
    strategies <- p$strategies
    pars <- strategies[[1]]$pars
    pars[[name]] <- value
    strategies[[1]]$pars <- pars
    p$strategies <- strategies
    ladder_assert_one_parameter(before, unlist(p$strategies[[1]]$pars), name)
    stand <- ladder_run(p)
    list(census = stand_census(stand),
         seed_height = ladder_as_patch(stand)$species[[1]]$new_node$height)
  }
  # The base value comes from the registered set, so a trait the strategy does not
  # carry as a parameter has no column here to referee.
  value <- unlist(build()$strategies[[1]]$pars)[[name]]
  h <- abs(value) * rel
  a <- census_at(value + h)
  b <- census_at(value - h)
  list(gradient = (a$census - b$census) / (2 * h),
       # The channel the differentiated paths impose to zero, measured rather than
       # assumed: preparation re-runs, so the seed height moves and this is by how
       # much.
       seed_height_slope = (a$seed_height - b$seed_height) / (2 * h))
}

# The same rebuilding difference on a stand where the species COMPETE.
#
# ladder_run_difference builds one species, and a one-species stand cannot see a
# whole class of defect: with no competition every cohort sits close to the state
# the leaf's own coefficients were read at, so a row that is only wrong far from
# there reads as correct. Measured -- the root-carbon rows agreed with a
# single-species reference to 1e-06 for both trait sets while the slow species'
# root allometry sat at 0.80 of this one.
#
# `species` picks which of the two the parameter is moved in, and the assertion
# that exactly one carried parameter moved is per species, so a hyperparameter
# reaching both is refused here as it is there.
ladder_run_difference_pair <- function(name, species = 2L,
                                       steps = c(1e-6, 1e-5, 1e-4, 1e-3),
                                       times = list(c(0, 0.63), c(0, 0.41)),
                                       regime = NULL) {
  # A regime is drivers AND, where it names one, a stand-wide extinction
  # coefficient, so it has to reach both halves of the build. Absent, this is the
  # default environment and the species' own traits.
  build <- function() {
    p <- ladder_parameters(c("fast", "slow"), k_I = regime$k_I)
    p$node_schedule_times <- times
    p
  }
  env_of <- function() {
    ladder_environment(regime$rain, if (is.null(regime$amplitude)) 0
                                   else regime$amplitude)
  }
  # ⚠️ BOTH SIDES CENSUS ON ONE TIME GRID, AND WITHOUT IT THIS MEASURES THE
  # SOLVER RATHER THAN THE MODEL. Moving a parameter by one part in a million
  # changes which steps the error estimator accepts, and the census then lands a
  # fixed distance away however small the move was -- so the quotient grows as
  # 1/h and every reading under that floor carries a derivative's units and none
  # of its content. DO NOT drop the pin to save the base run: unpinned, 1.theta's
  # leaf_area on drought reads 1.9e+05, -3.4e+04, -1196, -1043 over the four
  # steps, against -815.2, -815.7, -816.6, -744.4 pinned and -815.7 from the
  # sweep -- and no column on that stand did better, the census difference
  # tracking the step at a slope of 0.65 at best over its 86 readable ones where
  # a resolved difference gives 1.
  #
  # The TIMES only, not the step sizes with them. Pinning both makes each side
  # replay a step the other's state cannot take, which refuses outright on the
  # shaded and clamped stands; the grid is what the census accumulates on and is
  # all that has to agree.
  base <- ladder_run(build(), env = env_of())
  ode_times <- base$ode_times
  census_at <- function(value) {
    p <- build()
    strategies <- p$strategies
    before <- unlist(strategies[[species]]$pars)
    pars <- strategies[[species]]$pars
    pars[[name]] <- value
    strategies[[species]]$pars <- pars
    p$strategies <- strategies
    ladder_assert_one_parameter(before, unlist(p$strategies[[species]]$pars), name)
    p$ode_times <- ode_times
    stand_census(ladder_run(p, env = env_of()))
  }
  value <- unlist(build()$strategies[[species]]$pars)[[name]]
  at <- function(r) {
    h <- abs(value) * r
    (census_at(value + h) - census_at(value - h)) / (2 * h)
  }
  # ⚠️ The step-stability guard is not optional here, and skipping it once
  # produced a reading that reversed a decision. A column whose sensitivity is
  # an order below its neighbour's -- which competition produces, because the
  # suppressed species contributes less -- sits closer to this difference's own
  # floor, and the answer then moves with the step instead of holding. Measured
  # on the slow species' root allometry: 0.275, 0.131, 0.122, 0.126 across four
  # steps, so a single reading at 1e-05 is off by a factor of two from the
  # converged one.
  g <- vapply(steps, at, numeric(length(stand_census(base))))
  # The most-agreeing adjacent pair brackets where truncation crosses round-off;
  # its coarser member carries the less round-off of the two. One step for the
  # whole column, so a metric is not read at a different place from its neighbour.
  gaps <- abs(g[, -1L, drop = FALSE] - g[, -ncol(g), drop = FALSE])
  best <- which.min(apply(gaps, 2, max))
  # ⚠️ THE SPREAD IS ABSOLUTE, IN THE METRIC'S OWN UNITS, AND THE ONE CONSUMER
  # DIVIDES IT BY THE SAME SCALE IT DIVIDES ITS RESIDUAL BY. Normalising here too
  # is the second of two divisions and each hides the other: against the column's
  # own value a row at the difference's floor reads infinitely unconverged, and
  # against the largest of the four readings -- which on a series falling like
  # 1/h is its noisiest -- 79 of drought's 86 columns kept a reading further from
  # its neighbour than from zero and reported 0.03% to 34% for it.
  list(gradient = g[, best + 1L],
       column = paste0(species, ".", name),
       step = steps[[best + 1L]],
       values = g,
       spread = gaps[, best])
}

ladder_shortlist <- function() {
  # Two classes, and the second was unrefereed until a defect in it was found from
  # outside this ladder.
  #
  # Reduction-borne: each reaches a census through a reduction as well as through a
  # plant, so each can read zero for a structural reason rather than an ecological
  # one.
  #
  # Birth-size: each reaches a census through the seed's own size, so each carries
  # the one channel whose row is written where a newborn's state is. They move
  # together as a class -- a row missing there moved all eight while every other
  # channel held -- so the class is the unit, not the two of it that happened to be
  # listed here.
  unique(c("k_I", "eta", "a_l1", "a_l2", "a_st3", "recruitment_decay",
           ladder_birth_size_parameters()))
}

# The parameters that reach the census through birth size, where both
# differentiated paths impose the seed height's derivative to zero. They are the
# ones a whole-run difference is for, and they are named here rather than derived
# because the imposition is a declaration and not a property of the arithmetic:
# the seed height solves mass_live(h) = omega, and mass_live sums leaf, sapwood,
# bark and root mass.
#
# Reaching all eight needs a difference that perturbs a registered parameter
# rather than a trait, because only three of them are traits the fixture carries.
ladder_birth_size_parameters <- function() {
  c("lma", "rho", "omega", "theta", "a_l1", "a_l2", "a_r1", "a_b1")
}

# Which registered parameters a trait moves through the hyperparameter function,
# read by differencing that function rather than by reading it. Every entry
# besides the trait's own is a channel the sweep reports in its own column, and a
# difference taken on the trait vector would fold them all into one.
ladder_trait_fanout <- function(traits, name, rel = 1e-6) {
  derived <- function(tr) {
    unlist(as.data.frame(TF24_hyperpar(
      trait_matrix(unname(tr), names(tr)), TF24_Strategy(), filter = FALSE))[1, ])
  }
  h <- abs(traits[[name]]) * rel
  up <- traits; up[[name]] <- traits[[name]] + h
  dn <- traits; dn[[name]] <- traits[[name]] - h
  slope <- (derived(up) - derived(dn)) / (2 * h)
  slope[abs(slope) > 0 & names(slope) != name]
}

# Whether that difference is in its own domain on this fixture, which is not
# assumed. A re-run difference is unusable at production -- a relative step of two
# parts in ten million in leaf mass per area moves a mature stand between alive and
# identically zero -- so the check on the check is that the answer holds its figures
# across steps. Where it does not, the run is invalid rather than failing, as a
# violated regime assertion is.
#
# Four steps over three orders, and the value taken where they agree rather than at
# the finest. Truncation falls with the step and round-off rises, so the two ends
# fail in opposite directions and a series that reports the finest reports the
# round-off end: on a three-year run this reference reads 1.026 of its own plateau
# at a relative step of 1e-6. Three steps over two orders cannot say which end it is
# holding, and a spread taken over the whole series then reads round-off as drift.
#
# `plateau_at` is which end the agreement is at, and it belongs beside any figure
# taken from here.
ladder_run_difference_stable <- function(traits, name, ...,
                                         steps = c(1e-6, 1e-5, 1e-4, 1e-3)) {
  got <- lapply(steps, function(r) ladder_run_difference(traits, name, rel = r, ...))
  g <- vapply(got, function(x) x$gradient, numeric(length(got[[1]]$gradient)))
  scale <- max(abs(g))
  n <- length(steps)
  adjacent <- if (scale > 0) {
    vapply(seq_len(n - 1L),
           function(i) max(abs(g[, i] - g[, i + 1L])) / scale, numeric(1))
  } else rep(0, n - 1L)
  # The most-agreeing adjacent pair brackets the crossing of the two error terms;
  # take its coarser member, which carries the less round-off of the two.
  at <- which.min(adjacent)
  list(gradient = got[[at + 1L]]$gradient,
       seed_height_slope = got[[at + 1L]]$seed_height_slope,
       step = steps[[at + 1L]],
       steps = steps,
       values = g,
       adjacent = adjacent,
       plateau_at = if (at == 1L) "the small-step end"
                    else if (at == n - 1L) "the large-step end" else "the middle",
       spread = adjacent[[at]],
       spread_all = if (scale > 0) {
         max(apply(g, 1, function(r) diff(range(r)))) / scale
       } else 0)
}

# A unit direction in trait space, named as the gradient's columns are.
ladder_trait_direction <- function(columns, name) {
  at <- match(name, columns)
  stopifnot(!is.na(at))
  replace(numeric(length(columns)), at, 1)
}

# How far the replay's own census lands from the run's. This is the trajectory
# checks' floor: the reference replays step sizes rather than times, so its state
# at the end is the run's state only up to the arithmetic of re-stepping, and
# every derivative it reports is a derivative at the state it actually reached.
ladder_replay_floor <- function(stand, value) {
  reached <- ladder_trajectory_tangent(
    stand, numeric(length(census_trait_names_tf24(stand))))$value
  max(abs(reached - value) / abs(value))
}

# What the trajectory sweep and the trajectory tangent agree to. Every column is
# held here; there is no second, looser bound, and the history of the one there was
# is the argument for not having it.
#
# Measured worst 6.0e-05 over all 88 columns of a stand carrying one cohort per
# species. On the introductions fixture the largest is now recruitment_decay at
# 4.4e-05, with the two constants of the inflow condition next to it -- they sat at
# 1.9e-02 and 2.0e-03 while the boundary node's density adjoint was accumulated and
# never transposed through the condition that sets it.
#
# The two allometric constants were held to a separate bound of 4e-03 while a species
# carrying more than one cohort disagreed by 3.4e-03. That was a real defect and it is
# closed: the seed's dependent auxiliaries were not re-derived where its recorded height
# is written, so the leaf area every rate at birth size is scaled by carried no
# derivative. They now read 2.2e-09 and 1.1e-08 on the same fixture.
#
# Two things that bound cost, and they are why one bound is safer than two. It was the
# eight-parameter birth-size CLASS that was wrong, not the pair -- the other six were
# compared nowhere, so the bound named the symptom's two loudest columns and left the
# rest unrefereed. And a bound widened around a live disagreement cannot fail when that
# disagreement grows.
ladder_trajectory_agreement <- function() 3e-04

# What the transpose and the two references agree to on the soil rows, which is
# the tightest the leaf's supplied rows allow: a difference of the block re-solves
# the collar and the transpose records, so the gap is the rank-two factorisation's
# own fit step rather than a missing channel.
#
# It was 1.06e-05 while the uptake reduction's transpose scattered onto the
# introduced nodes only, leaving the boundary node's own draw -- the distribution's
# lower grid point -- reaching no soil row. Both reductions now carry that row and
# the residual is 5.38e-08 on the four-node patch, 5.07e-08 crossed and 4.67e-08
# on the one-cohort patch.
ladder_soil_row_agreement <- function() 8e-08

# A parameter the boundary cannot answer for. A refusal, never a number, and an
# unknown parameter must refuse by name rather than return anything at all.
ladder_refused_by_name <- function() {
  # The price of a cost curve this strategy does not run. TF24_floor is TF24 plus
  # a floor on the marginal value of water, and its lambda_o is that floor: the
  # leaf carries the parameter because the curve exists, and no equation TF24
  # evaluates reads it. So the column is exactly zero, and it is zero for a reason
  # that is a property of WHICH CURVE is seated rather than of the model's
  # sensitivity -- which is why it is refused by name rather than declared as an
  # insensitivity.
  #
  # `p_50` was here and is gone: it is `stem_P50` now, it is an ordinary trait of
  # the vulnerability curve, and its column is live.
  c("TF24_floor_lambda_o")
}

# The band a registered parameter that reaches no equation comes back in. It is
# round-off rather than exact zero, so the check is two-sided and not a test
# against zero.
ladder_roundoff_band <- function() c(1e-22, 1e-16)
