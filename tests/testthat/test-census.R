# The census reduction and the R-facing gradient entry point.

# TF24's allometry, written out from its equations rather than called through
# the strategy, so the value gate does not check the reduction against itself.
tf24_allometry_r <- function(pars, eta) {
  eta_c <- 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta)
  function(height, area_heartwood, mass_heartwood) {
    area_leaf <- (height / pars$a_l1)^(1 / pars$a_l2)
    area_sapwood <- area_leaf * pars$theta
    area_bark <- pars$a_b1 * area_leaf * pars$theta
    mass_leaf <- area_leaf * pars$lma
    mass_sapwood <- area_sapwood * height * eta_c * pars$rho
    mass_bark <- area_bark * height * eta_c * pars$rho
    list(leaf_area = area_leaf,
         mass_above_ground = mass_leaf + mass_bark + mass_sapwood +
           mass_heartwood,
         area_stem = area_bark + area_sapwood + area_heartwood)
  }
}

# Every cohort's state, boundary node last, from the ODE state each node writes.
# The density is the one each node contributes, which on the birth-date
# coordinate includes, for the newest node and the boundary node, the
# establishment of the interval still open; the birth date is not ODE state.
# Both are read off the nodes themselves.
species_state_r <- function(species) {
  nodes <- c(species$nodes, list(species$new_node))
  names <- species$new_node$ode_names
  s <- vapply(nodes, function(n) n$ode_state, numeric(length(names)))
  rownames(s) <- names
  data.frame(height = s["height", ],
             area_heartwood = s["area_heartwood", ],
             mass_heartwood = s["mass_heartwood", ],
             log_density = vapply(nodes, function(n) n$log_density, 0),
             birth_date = c(vapply(species$nodes,
                                   function(n) n$introduction_time, 0),
                            species$new_node$introduction_time))
}

trapezium_r <- function(x, y) {
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1))) / 2
}

# The census in R. A census is a quadrature of a density, so the grid is the
# coordinate the density is carried in: the birth date, which ascends as the
# nodes are stored and puts the boundary node last, or the height, which ascends
# the other way and puts it first.
census_r <- function(species, pars, eta, include_boundary = TRUE,
                     birth_date = TRUE) {
  st <- species_state_r(species)
  n <- nrow(st)
  order <- if (birth_date) {
    if (include_boundary) seq_len(n) else seq_len(n - 1)
  } else {
    if (include_boundary) c(n, seq(n - 1, 1)) else seq(n - 1, 1)
  }
  st <- st[order, , drop = FALSE]
  psi <- tf24_allometry_r(pars, eta)(st$height, st$area_heartwood,
                                     st$mass_heartwood)
  density <- exp(st$log_density)
  x <- if (birth_date) st$birth_date else st$height
  vapply(psi, function(p) trapezium_r(x, density * p), numeric(1))
}

# One stand for every block below, on the birth-date coordinate the reverse
# sweep transposes, and WRITTEN into the initial state rather than grown into.
# Three cohorts at mutually non-commensurate birth dates, three scheduled
# introductions and a patch lifetime of 2: six nodes and a live boundary node.
#
# ⚠️ THE LIVE BOTTOM IS THE PREMISE EVERY CHECK HERE RESTS ON, and it is what a
# grown distribution is expensive to reach. Where the boundary node's density has
# underflowed to zero so has its neighbour's, the closing trapezium contributes
# exactly nothing, and the reduction gate and the adjoint gates below go quiet
# rather than failing. Measured: the bottom two densities are 9.96e-01 and
# 9.92e-01, where the same recipe at a lifetime of 3 gives zero and a distribution
# grown to a lifetime of 5 gives zero.
#
# ⚠️ SIX NODES RATHER THAN THE SEVENTY-EIGHT A DEFAULT SCHEDULE FILLS IN IS WHAT
# MAKES THE BOUNDARY NODE READABLE. It is the reduction's closing grid point and
# the one part of the adjoint no state column reports, so the birth-date-grid
# test reads it as the gap between two references -- and that gap is 2.5e-07 of
# the main term here against 6.1e-09 under a default schedule, where the
# reference's own floor is 9.5e-11.
#
# What it does not carry is the claim that a real trajectory REACHES such a
# state. Nothing here claims that; the premise is a guard on vacuity.
census_stand <- local({
  built <- NULL
  function() {
    if (is.null(built)) {
      ctrl <- Control(node_density_coordinate = "birth_date")
      p <- testing_parameters("TF24")
      p$max_patch_lifetime <- 2
      p <- add_strategies(p, trait_matrix(0.0825, "lma"))
      state <- make_initial_state(
        p, heights = list(c(6.2, 3.1, 1.37)),
        log_densities = list(c(-1.3, -0.7, -0.21)),
        env = Environment("TF24"), ctrl = ctrl)
      state$node_times <- list(c(-0.41, -0.23, -0.07))
      p <- set_initial_state(p, state)
      p$node_schedule_times <- list(c(0, 0.73, 1.41))
      solver <- DeterministicPatchRunner("TF24", "TF24_Env")(p, Environment("TF24"), empty_events(),
                                     ctrl)
      solver$run()
      built <<- solver
    }
    built
  }
})

# The leaf-area census as a function of one entry of the ODE state, for
# differencing against the adjoint. `solver$patch` is a copy, so the perturbation
# cannot reach the stand every other block shares.
#
# ⚠️ THE BOUNDARY NODE HAS TO MOVE WITH THE STATE, and a reference that holds it
# fixed differentiates a different function. It is not ODE state:
# census_state_and_trait_rows rebuilds it through set_state_and_boundary, from a
# light field every cohort's height and density enters by the same product the
# census integrates. Holding it fixed leaves a gap that grows with the stand --
# 2.5e-07 here, 8.6e-05 on one grown to a lifetime of 12 -- and the
# birth-date-grid test below measures it rather than absorbing it in a tolerance.
census_of_state <- function(solver, column, birth_date = TRUE) {
  strategy <- solver$parameters$strategies[[1]]
  patch <- solver$patch
  y <- patch$ode_state
  time <- patch$time
  function(value) {
    y[[column]] <- value
    patch$set_ode_state(y, time)
    # A rate evaluation owns the field build and the inflow condition, so this
    # leaves the boundary node where set_state_and_boundary would.
    invisible(patch$ode_rates)
    census_r(patch$species[[1]], strategy$pars,
             strategy$pars$eta, birth_date = birth_date)[["leaf_area"]]
  }
}

test_that("the census value matches an independent R reduction", {
  solver <- census_stand()
  species <- solver$patch$species[[1]]
  strategy <- solver$parameters$strategies[[1]]
  expected <- census_r(species, strategy$pars, strategy$pars$eta)
  got <- stand_census(solver)
  expect_equal(names(got), c("leaf_area", "mass_above_ground", "area_stem"))
  expect_equal(unname(got), unname(expected[names(got)]), tolerance = 1e-12)
})

test_that("the census reduction includes the boundary node", {
  solver <- census_stand()
  species <- solver$patch$species[[1]]
  strategy <- solver$parameters$strategies[[1]]
  st <- species_state_r(species)
  # The premise: the closing interval has live density at both ends.
  expect_gt(exp(st$log_density[nrow(st)]), 0)
  expect_gt(exp(st$log_density[nrow(st) - 1]), 0)
  with_boundary <- census_r(species, strategy$pars, strategy$pars$eta, TRUE)
  without <- census_r(species, strategy$pars, strategy$pars$eta, FALSE)
  # A reduction that starts at the smallest cohort drops the interval down to
  # the boundary node. Measured, the closing interval is 3.5e-04 of the answer.
  expect_gt(min(abs(with_boundary - without) / abs(with_boundary)), 1e-6)
  expect_equal(unname(stand_census(solver)), unname(with_boundary),
               tolerance = 1e-12)
})

test_that("the state adjoint is the census's derivative on the birth-date grid", {
  # A height is state and a birth date is not, so which one the grid is built
  # from decides whether a cohort's height moves the quadrature as well as the
  # integrand. On this coordinate it does not.
  #
  # Every node and both state families, against a Richardson difference of the
  # census itself rather than of a reduction written out here: the grid claim is
  # about the abscissa, and a reference that re-derives the integrand as well is
  # answering two questions with one number.
  #
  # NOT gated on a refusal. This fixture is written rather than reached, so the
  # adjoint being refused here is the model having moved under a fixture that did
  # not -- which is the finding, and a skip would report it as green.
  solver <- census_stand()
  species <- solver$patch$species[[1]]
  strategy <- solver$parameters$strategies[[1]]
  st <- species_state_r(species)
  adjoint <- stand_census_state_adjoint(solver)
  names_i <- species$nodes[[1]]$ode_names
  stride <- length(names_i)
  n_node <- length(species$nodes)

  # The reference is taken on a copy of the patch, so the copy must start where
  # the stand is: a round trip through the state setter that moved anything
  # would put every difference below at a different point from the adjoint.
  expect_equal(census_of_state(solver, 1L)(solver$patch$ode_state[[1]]),
               unname(stand_census(solver)[["leaf_area"]]), tolerance = 1e-14)

  psi_at <- tf24_allometry_r(strategy$pars, strategy$pars$eta)
  b <- st$birth_date
  w <- vapply(seq_along(b), function(i) {
    lo <- if (i > 1) (b[i] - b[i - 1]) / 2 else 0
    hi <- if (i < length(b)) (b[i + 1] - b[i]) / 2 else 0
    lo + hi
  }, 0)
  density <- exp(st$log_density)
  area_leaf <- psi_at(st$height, st$area_heartwood, st$mass_heartwood)$leaf_area
  d_area_leaf <- (1 / strategy$pars$a_l2) *
    (st$height / strategy$pars$a_l1)^(1 / strategy$pars$a_l2 - 1) /
    strategy$pars$a_l1

  got <- boundary_held <- over_birth_date <- over_height <- numeric(0)
  for (slot in c("height", "log_density")) {
    for (k in seq_len(n_node)) {
      col <- (k - 1) * stride + match(slot, names_i)
      at <- solver$patch$ode_state[[col]]
      step <- 1e-2 * max(abs(at), 1)
      over_birth_date <- c(over_birth_date,
        test_gradient_richardson(census_of_state(solver, col), at, step, 4))
      # Rows run in storage order with the boundary node last, which ascends in
      # birth date and descends in height; negating the height integral is what
      # takes it up its own axis instead.
      over_height <- c(over_height,
        -test_gradient_richardson(census_of_state(solver, col, FALSE), at, step, 4))
      got <- c(got, adjoint["leaf_area", col])
      # The same derivative with the boundary node pinned, which is what a
      # reduction written over a fixed state table gives.
      boundary_held <- c(boundary_held, w[k] * density[k] *
        if (slot == "height") d_area_leaf[k] else area_leaf[k])
    }
  }

  # The adjoint is the derivative on the coordinate the density is carried in.
  # Measured, the worst node reads 9.5e-11 and most read below 1e-13.
  expect_lt(max(abs(got - over_birth_date) / abs(over_birth_date)), 1e-8)
  # And it is not the derivative on the other grid. Stated per node rather than
  # as one non-vacuity check: the two grids are 1.4x apart at the closest node
  # and 54x at the furthest, so an adjoint built on the wrong one is nowhere a small
  # error.
  expect_gt(min(abs(got - over_height) / abs(over_birth_date)), 0.25)
  # The boundary node does not respond to these columns. Its density is the
  # birth rate times the establishment the open interval has integrated, which
  # is ODE state of its own; the field reaches establishment through that
  # state's rate, not through the census. So the derivative with the boundary
  # node pinned is the whole answer here, and a census that read establishment
  # off the field at the instant would show up as a gap.
  expect_lt(max(abs(got - boundary_held) / abs(got)), 1e-10)
})

test_that("the state adjoint reaches every state a metric reads", {
  # An exact zero in this design is the signature of a missing accumulator and
  # never of true insensitivity, so a state a metric demonstrably reads must have
  # an adjoint.
  #
  # `area_stem` sums bark, sapwood and heartwood AREA, and `mass_above_ground`
  # sums the three masses and heartwood MASS, so both read a heartwood state
  # directly and linearly: the adjoint is the quadrature weight times the density,
  # with no allometry in between, which is what makes this checkable by hand.
  # The R-reduction test shows the census VALUE carries them, so the reduction
  # is right and it is the recording of it that would not be.
  #
  # Exact rather than approximate, and the reason is the finding the
  # birth-date-grid test rests on: a heartwood state is not in the light field,
  # so these two columns carry no boundary-node channel and nothing else stands
  # between the weight and the answer.
  solver <- census_stand()
  adjoint <- stand_census_state_adjoint(solver)
  species <- solver$patch$species[[1]]
  names_i <- species$nodes[[1]]$ode_names
  stride <- length(names_i)

  # Birth-date trapezium weights, boundary node last, as the census integrates.
  b <- c(vapply(species$nodes, function(n) n$introduction_time, 0),
         species$new_node$introduction_time)
  dens <- c(vapply(species$nodes, function(n) exp(n$log_density), 0),
            exp(species$new_node$log_density))
  w <- vapply(seq_along(b), function(i) {
    lo <- if (i > 1) (b[i] - b[i - 1]) / 2 else 0
    hi <- if (i < length(b)) (b[i + 1] - b[i]) / 2 else 0
    lo + hi
  }, 0)
  expected <- (w * dens)[seq_along(species$nodes)]

  expect_true(all(expected > 0))
  for (pair in list(c("area_stem", "area_heartwood"),
                    c("mass_above_ground", "mass_heartwood"))) {
    cols <- (seq_along(species$nodes) - 1) * stride + match(pair[2], names_i)
    expect_equal(unname(adjoint[pair[1], cols]), expected,
                 tolerance = 1e-12, info = paste(pair[1], "reads", pair[2]))
  }

  # Non-vacuity: the states that DO have an adjoint, so a wholesale failure of the
  # recording would not pass this by looking like the defect above.
  h_cols <- (seq_along(species$nodes) - 1) * stride + match("height", names_i)
  expect_true(any(adjoint["leaf_area", h_cols] != 0))
})

test_that("stand_gradient_compare() refuses to compare across two Controls", {
  a <- list(gradient = matrix(1, 1, 1, dimnames = list("leaf_area", "lma")),
            control = c(GSS_tol_abs = 1e-3, ci_abs_tol = 1e-3,
                        node_gradient_eps = 1e-6, schedule_eps = 1e-3,
                        gradient_curvature_floor = 1e-3))
  b <- a
  b$control[["schedule_eps"]] <- 1e-4
  expect_error(stand_gradient_compare(a, b), "schedule_eps")
  expect_equal(stand_gradient_compare(a, a),
               matrix(0, 1, 1, dimnames = list("leaf_area", "lma")))
})

test_that("no census metric has an all-zero state sensitivity", {
  # A structural claim about the adjoint rather than one about the shape of the
  # distribution, so it reads the stand every other block uses.
  solver <- census_stand()
  # The census's rows; the fitness rows read the fecundity state.
  adjoint <- stand_census_state_adjoint(solver)[names(stand_census(solver)), ,
                                          drop = FALSE]
  expect_equal(nrow(adjoint), 3L)
  # Every metric is built from height, and every cohort's log density multiplies
  # it, so both state families must move all three metrics.
  stride <- length(solver$patch$species[[1]]$new_node$ode_names)
  n_node <- solver$patch$species[[1]]$size
  height_cols <- (seq_len(n_node) - 1) * stride + 1
  density_cols <- (seq_len(n_node) - 1) * stride + stride
  for (m in rownames(adjoint)) {
    expect_true(any(adjoint[m, height_cols] != 0))
    expect_true(any(adjoint[m, density_cols] != 0))
  }
})

test_that("the Control a gradient is taken at is the entries that move it", {
  # Five move the TRAJECTORY the gradient is taken along. The sixth moves no
  # forward number at all and still decides which rows exist, by refusing a
  # collar response the profit curvature is too small to support -- so two
  # gradients taken at different floors are gradients of different functions for a
  # different reason, and both reasons belong in the same comparison.
  #
  # ⚠️ ASKED BY NAME, WHICH IS WHAT MAKES THIS ABLE TO FAIL. A list of names kept
  # in R, with values compared POSITIONALLY against a list built in the same
  # order, is self-consistent either way. Asked by name, a value wired to the
  # wrong Control field FAILS, and so does a name changed in
  # deterministic_patch_runner.h.
  #
  # ⚠️ What it still cannot see is a swap among GSS_tol_abs, ci_abs_tol and
  # gradient_curvature_floor, because all three are 1e-3 at the defaults and a swap
  # of equal numbers has nothing to read. That is fine HERE and only here: each
  # travels from deterministic_patch_runner.h as name-with-value, so there is no
  # second ordering left to
  # disagree with and nothing can produce that swap any more. Do not answer it
  # with a test.
  solver <- census_stand()
  got <- gradient_control(solver)
  expect_equal(names(got),
               c("GSS_tol_abs", "ci_abs_tol", "node_gradient_eps",
                 "schedule_eps", "gradient_curvature_floor",
                 "offspring_height_tol"))
  expect_equal(got, unlist(unclass(Control())[names(got)]))
  # A newborn's height tolerance moves every model's trajectory, so two
  # gradients taken at different ones are not compared.
  a <- list(control = got, gradient = matrix(0))
  b <- a
  b$control[["offspring_height_tol"]] <- 1e-8
  expect_error(stand_gradient_compare(a, b), "offspring_height_tol")
})

test_that("the trait gradient entry point is reachable", {
  # The symbol exists, stand_gradient reaches it, and on this stand it answers:
  # a refusal here is the model having moved under a written fixture, which is
  # the finding rather than a state to tolerate.
  solver <- census_stand()
  expect_true(is.function(census_trait_gradient))
  # A column is named for its species as well as its parameter: a bare name would
  # resolve to species one's column silently on a multi-species stand.
  got <- stand_gradient(solver, traits = "1.lma", hyperpar = NULL)
  expect_equal(rownames(got$gradient), unname(census_metric_names(solver, "TF24")))
  expect_equal(colnames(got$gradient), "1.lma")
  expect_false(any(stand_gradient_refused(got)))
  expect_true(all(is.finite(got$gradient)))
  # And the bare name refuses, naming the convention rather than only failing.
  expect_error(stand_gradient(solver, traits = "lma", hyperpar = NULL), "species index")
})

test_that("the gradient refuses, by name, what the strategy cannot differentiate", {
  short <- function(edit = identity, ctrl = Control()) {
    p <- testing_parameters("TF24")
    p$max_patch_lifetime <- 0.3
    p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
    p$node_schedule_times <- list(c(0, 0.1, 0.2))
    p <- edit(p)
    run(p, Environment("TF24"), ctrl, record = "states")$solver
  }
  # Each is a property of the strategy alone, so it is refused before the run
  # is recorded, naming the entry point and the species.
  # Deep-crown shading's crown means pass through the leaf at double.
  solver <- short(ctrl = Control(shading_model = "deep-crown"))
  expect_error(stand_gradient(solver, metrics = "leaf_area", hyperpar = NULL),
               "species 1: shading_model 'deep-crown' is not differentiable")
  # Sapwood acclimation's sensitivity is a difference through the leaf, which
  # carries no rows.
  solver <- short(function(p) {
    p$strategies[[1]]$pars$a_sw <- 0.1
    p
  })
  expect_error(stand_gradient(solver, hyperpar = NULL),
               "census_trait_gradient: species 1: sapwood acclimation \\(a_sw > 0\\) is not differentiable")
  # The tangent reference differentiates the same rates, and is refused alike
  # rather than answering without the controller's sensitivity.
  direction <- rep(1, length(census_trait_names(solver, "TF24")))
  expect_error(ladder_trajectory_tangent_tf24(solver, direction),
               "census_trait_tangent: species 1: sapwood acclimation")
  # The census's state adjoint at one state reads no sapwood sensitivity, and
  # still answers.
  expect_true(all(is.finite(stand_census_state_adjoint(solver))))
})

test_that("TF24w has no stand gradient", {
  # Its private root-zone columns couple cohorts outside the recorded block, so
  # no census is declared for it, and no gradient entry point is compiled.
  p <- testing_parameters("TF24w")
  p$max_patch_lifetime <- 0.3
  p <- add_strategies(p, trait_matrix(0.0825, "lma"),
                      hyperpar = default_hyperpar("TF24"), birth_rate = 20)
  p$node_schedule_times <- list(c(0, 0.1, 0.2))
  solver <- run(p, Environment("TF24w"),
                control(node_density_coordinate = "birth_date"))$solver
  expect_error(stand_gradient(solver, hyperpar = NULL),
               "No census is declared for the TF24w model")
})

test_that("an empty request and an incomparable pair are refused", {
  solver <- census_stand()
  expect_error(stand_gradient(solver, metrics = character(0), hyperpar = NULL),
               "names no metric")
  a <- stand_gradient(solver, metrics = "leaf_area", traits = "1.lma",
                      hyperpar = NULL)
  b <- a
  b$control <- b$control[-1]
  expect_error(stand_gradient_compare(a, b), "different Control entries")
  expect_equal(stand_gradient_compare(a, a), a$gradient - a$gradient)
})
