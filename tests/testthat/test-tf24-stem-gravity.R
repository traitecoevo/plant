# The stem's gravitational head, rho*g*H*eta_c, gated by use_stem_gravity.
#
# Gravity is a separate term from the stem's resistance, not part of it: it is
# paid at zero flow, where a resistance costs nothing. So the test that matters
# is the predawn one -- a shut leaf on wet soil sits a head above the collar,
# and the head is the lift to the leaf-area-weighted mean leaf height. See
# TF24_Model::soil_potential_at_leaf_height for where the head enters and what
# it costs to treat it as one offset.

# rho*g as phylloptim carries it, MPa per metre: the soil layers' heads use the
# same constant, so the stem's is pinned to it here.
GRAVITY_HEAD <- 9.8e-3

# One leaf at `height` over one soil layer at volumetric water content `theta`,
# under `SWdown` W m^-2: the auxiliary readings, with the soil potential and the
# layer midpoint the leaf was handed.
gravity_probe <- function(strategy, height, theta = 0.25, SWdown = 900) {
  ind <- TF24_Individual(strategy)
  ind$set_state("height", height)
  env <- TF24_Environment()
  env$set_soil_number_of_depths(1)
  env$set_soil_water_state(theta)
  env$extrinsic_drivers_set_constant("SWdown", SWdown)
  ind$compute_rates(env)
  c(stats::setNames(ind$internals$auxs, ind$aux_names),
    psi_soil = env$psi_from_soil_moist(theta),
    z_mid = env$get_soil_mid_depths()[[1]])
}

with_gravity <- function(s = TF24_Model()) {
  s$pars$use_stem_gravity <- 1
  s
}

test_that("the stem's head is off by default", {
  expect_identical(TF24_Model()$pars$use_stem_gravity, 0)
})

test_that("a shut leaf sits the stem's head above the collar", {
  # In the dark the leaf has no reason to open, so all but a trace of flow
  # stops and every resistance drops out: the leaf potential is the soil's
  # plus the lift from the layer to the leaf, which is the predawn gradient.
  #
  # The trace is ~1e-13 of transpiration, the solve's own tolerance on a shut
  # leaf, and it leaves the collar ~2.5e-07 MPa from the soil's balance at
  # 20 m. The gap between leaf and collar is the head itself to ~5e-09, so it
  # is held tighter than the potentials it is the difference of.
  s <- with_gravity()
  eta_c <- 1 - 2 / (1 + s$pars$eta) + 1 / (1 + 2 * s$pars$eta)
  for (h in c(1, 5, 20)) {
    on <- gravity_probe(s, h, SWdown = 0)
    off <- gravity_probe(TF24_Model(), h, SWdown = 0)
    expect_lt(on[["transpiration"]], 1e-12)
    expect_lt(off[["transpiration"]], 1e-12)
    head <- GRAVITY_HEAD * h * eta_c
    balance <- on[["psi_soil"]] + GRAVITY_HEAD * on[["z_mid"]]
    expect_lt(abs(on[["opt_psi_leaf"]] - (balance + head)), 1e-6,
              label = paste("leaf at", h, "m"))
    # The reported collar is the real one, below the head: at the soil's
    # balance with the head on as with it off, and the off leaf sits at it.
    expect_lt(abs(on[["opt_root_psi"]] - balance), 1e-6)
    expect_lt(abs(off[["opt_psi_leaf"]] - off[["opt_root_psi"]]), 1e-8)
    expect_lt(abs(on[["opt_psi_leaf"]] - on[["opt_root_psi"]] - head), 1e-7,
              label = paste("head at", h, "m"))
  }
})

test_that("with water moving, the leaf pays the head on top of the friction", {
  s <- with_gravity()
  eta_c <- 1 - 2 / (1 + s$pars$eta) + 1 / (1 + 2 * s$pars$eta)
  for (h in c(2, 10)) {
    on <- gravity_probe(s, h)
    off <- gravity_probe(TF24_Model(), h)
    expect_gt(off[["transpiration"]], 0)
    # Drier leaf, less water, less carbon.
    expect_gt(on[["opt_psi_leaf"]], off[["opt_psi_leaf"]])
    expect_lt(on[["transpiration"]], off[["transpiration"]])
    expect_lt(on[["assimilation"]], off[["assimilation"]])
    # Leaf below collar by friction plus the head, and the friction is not
    # negative.
    expect_gt(on[["opt_psi_leaf"]] - on[["opt_root_psi"]],
              GRAVITY_HEAD * h * eta_c)
  }
})

test_that("the head costs a tall plant more than a short one", {
  # The head is linear in height while the friction, under conduit widening,
  # grows as H^(1 - beta): so the share of carbon gravity takes must rise with
  # height. A constant penalty would mean the head never reached the leaf.
  s <- with_gravity()
  cost <- function(h) {
    1 - gravity_probe(s, h)[["assimilation"]] /
      gravity_probe(TF24_Model(), h)[["assimilation"]]
  }
  expect_gt(cost(1), 0)
  expect_gt(cost(15), cost(1))
})

test_that("with the head on, a growth trait's column is the pinned-step difference", {
  # The head is formed from height at the strategy's scalar, so it carries
  # height's rows into the leaf's supply rows. Formed at double instead, every
  # column misses the head's share of the response through height: built that
  # way, all four comparisons below fail.
  p <- tf24_fitness_stand_parameters()
  p$strategies[[1]]$pars$use_stem_gravity <- 1
  solver <- tf24_fitness_stand_run(p)
  off <- tf24_fitness_stand_parameters()
  expect_false(isTRUE(all.equal(
    stand_fitness(tf24_fitness_stand_run(off, record = FALSE))[[1]],
    stand_fitness(solver)[[1]], tolerance = 1e-6)))

  metrics <- c("1.net_reproduction_ratio", "1.offspring_production")
  names <- c("hmat", "K_s")
  g <- stand_gradient(solver, metrics = metrics, traits = paste0("1.", names),
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))

  pinned <- pin_steps(solver)
  fitness_of <- function(q) stand_fitness(tf24_fitness_stand_run(q, record = FALSE))
  for (name in names) {
    fd <- central_difference(pinned, 1, name, NULL, fitness_of, rel = 1e-5)
    for (metric in metrics) {
      expect_lt(abs(g$gradient[metric, paste0("1.", name)] - fd[[metric]]),
                1e-5 * abs(fd[[metric]]), label = paste(metric, name))
    }
  }
})
