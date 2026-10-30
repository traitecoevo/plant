# TF24w: TF24 with a private root-zone soil column per cohort (#681). The model
# is in inst/include/plant/models/tf24w_strategy.h and the coupling passes in
# Patch::private_water_aggregates / private_water_exchange.

# A short pulsed stand, the paper's settings: one lma, annual cohorts, rain as
# Poisson pulses from a fixed seed, the continuous driver off, the birth-date
# coordinate.
tf24w_stand <- function(type, mixing = 0, years = 6, birth_rate = 5,
                        rainfall_pulses = TRUE, absorb = 0, seed = 1,
                        schedule = seq(0, years - 1, by = 1), lma = 0.0825) {
  p <- scm_base_parameters(type)
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- years
  p <- add_strategies(p, trait_matrix(lma, "lma"), hyperpar = hyperpar(type),
                      birth_rate = birth_rate)
  p$node_schedule_times <- list(schedule)
  env <- Environment(type)
  env$extrinsic_drivers_set_constant("rainfall", 0)
  env$water_mixing_rate <- mixing
  if (rainfall_pulses) {
    set.seed(seed)
    s <- pulse_series_poisson(0, years, rate = 20, mean_amount = 0.014)
    ev <- events(node_introductions(p), rainfall_pulse(s$time, s$amount))
  } else {
    ev <- events(node_introductions(p))
  }
  ctrl <- control(node_density_coordinate = "birth_date", ode_step_size_max = 0.05,
                  absorb_epsilon = absorb)
  run_scm(p, env, ctrl, events = ev, run = FALSE)
}

test_that("TF24w declares its states, readings and version", {
  s <- TF24w_Strategy()
  nm <- TF24w_Individual(s)$ode_names
  expect_identical(nm, c(TF24_Individual()$ode_names, paste0("water_", 1:5)))
  expect_identical(model_id("TF24w"), paste0("TF24w@v", model_version("TF24"), ".1"))
  expect_true(all(c("water_infiltration", "water_drainage", "crown_area") %in%
                    TF24w_Individual(s)$aux_names))
  # The crown-LAI parameter is TF24's, inert there, and reaches the hyperpar.
  expect_true("L_crown" %in% names(s$pars))
  expect_true("L_crown" %in% names(TF24_Strategy()$pars))
  env <- Environment("TF24w")
  expect_identical(env$water_mixing_rate, 0)
})

test_that("TF24w refuses the height coordinate", {
  p <- scm_base_parameters("TF24w")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = TF24w_hyperpar)
  expect_error(run_scm(p, Environment("TF24w"),
                       control(node_density_coordinate = "height"), run = FALSE),
               "birth-date")
})

# Each node's water column, as a matrix of layers by node, and the environment's.
tf24w_columns <- function(scm) {
  sp <- scm$patch$species[[1]]
  nodes <- sp$nodes
  nm <- nodes[[1]]$ode_names
  first <- which(nm == "water_1")
  list(water = sapply(nodes, function(n) n$ode_state[first:(first + 4)]),
       log_density = sapply(nodes, function(n) n$ode_state[which(nm == "log_density")]),
       height = sapply(nodes, function(n) n$ode_state[which(nm == "height")]),
       bare = scm$patch$environment$get_soil_water_state())
}

test_that("every column holds a water content inside the soil's range", {
  scm <- tf24w_stand("TF24w", years = 3)
  scm$run_to(2.5)
  cols <- tf24w_columns(scm)
  expect_true(all(cols$water >= 0))
  expect_true(all(cols$water <= 0.428))
  expect_equal(nrow(cols$water), 5)
})

test_that("with no mixing a cohort dries its own column and leaves the bare ground wet", {
  # A single cohort at a vanishing density: on TF24 it could never dry the patch
  # soil, which is the zero-density problem #681 exists to remove. Here its
  # column is its own, however rare it is. No rain, so the bare column only
  # drains.
  dry <- function(mixing) {
    scm <- tf24w_stand("TF24w", mixing = mixing, years = 2, birth_rate = 1e-6,
                       rainfall_pulses = FALSE, schedule = 0)
    scm$run_to(2)
    cols <- tf24w_columns(scm)
    list(own = cols$water[1, 1], bare = cols$bare[1])
  }
  private <- dry(0)
  expect_lt(private$own, 0.75 * private$bare)
  pooled <- dry(200)
  expect_equal(pooled$own, pooled$bare, tolerance = 0.1)
  expect_lt(abs(pooled$own - pooled$bare), 0.2 * abs(private$own - private$bare))
})

test_that("the pooled limit reproduces TF24", {
  # At a large mixing rate every column holds the patch mean, rain and drainage
  # per unit area are the same everywhere, and the summed uptake per m2 is
  # TF24's depletion: the two models should agree to the mixing lag.
  ref <- tf24w_stand("TF24", years = 6)
  ref$run()
  w <- tf24w_stand("TF24w", mixing = 200, years = 6)
  w$run()
  # scm$patch$environment is a copy, so it is read after the run.
  ew <- w$patch$environment
  er <- ref$patch$environment
  expect_equal(ew$get_soil_water_state(), er$get_soil_water_state(),
               tolerance = 0.03)
  expect_equal(w$net_reproduction_ratios, ref$net_reproduction_ratios,
               tolerance = 0.03)
  # And the accumulated budget terms, area-weighted over every column in TF24w,
  # agree with the one column of TF24.
  fw <- ew$get_soil_water_state_cumulative_flux()
  fr <- er$get_soil_water_state_cumulative_flux()
  expect_gt(fr[2], 0)
  expect_equal(fw, fr, tolerance = 0.03)
})

test_that("the resident's own copy, introduced as a rare second species, scores exactly 0", {
  p <- scm_base_parameters("TF24w")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 4
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.0825), "lma"),
                      hyperpar = TF24w_hyperpar, birth_rate = list(5, 5e-6))
  p$node_schedule_times <- list(0:3, 0:3)
  env <- Environment("TF24w")
  env$extrinsic_drivers_set_constant("rainfall", 0.3)
  env$water_mixing_rate <- 1
  scm <- run_scm(p, env, control(node_density_coordinate = "birth_date",
                                 ode_step_size_max = 0.05), run = FALSE)
  scm$run()
  R <- scm$net_reproduction_ratios
  expect_equal(R[[2]], R[[1]], tolerance = 1e-10)
})

# The patch's water, as the area-weighted total over every column, in metres:
# each node's share is its crown area (from its height through the leaf-area
# allometry and L_crown) times its density times its trapezium weight on the
# birth-date grid, the boundary node closing the grid but owning no ground;
# shares are scaled where crowns overlap; the bare ground is what crowns leave.
tf24w_storage <- function(scm, dz = 0.3) {
  sp <- scm$patch$species[[1]]
  s <- scm$patch$parameters$strategies[[1]]$pars
  crown <- function(h) (h / s$a_l1)^(1 / s$a_l2) / s$L_crown
  cols <- tf24w_columns(scm)
  x <- c(sp$node_times, scm$patch$time)
  n <- exp(cols$log_density)
  h <- cols$height
  water <- cols$water
  m <- length(x)
  w <- numeric(m - 1)
  for (k in seq_len(m - 1)) {
    left <- if (k > 1) x[k] - x[k - 1] else 0
    right <- x[k + 1] - x[k]
    w[k] <- (left + right) / 2
  }
  share <- w * n * crown(h)
  cover <- sum(share)
  share <- share / max(1, cover)
  bare <- max(0, 1 - cover)
  list(total = sum(share * colSums(water) * dz) + bare * sum(cols$bare * dz),
       cover = cover)
}

test_that("the water budget closes over births, growth, death and pulses", {
  # Storage change equals infiltration (pulses included) minus drainage minus
  # uptake, all area-weighted over the columns: the land-transfer terms are what
  # make this hold while shares move.
  scm <- tf24w_stand("TF24w", mixing = 1, years = 12, birth_rate = 20)
  # Every column starts at the environment's content. scm$patch$environment is
  # a copy, so the fluxes are read from a fresh handle after the run.
  start <- sum(scm$patch$environment$get_soil_water_state() * 0.3)
  scm$run()
  flux <- scm$patch$environment$get_soil_water_state_cumulative_flux()
  st <- tf24w_storage(scm)
  # sum_rainfall, sum_infiltration, sum_drainage, sum_resource_depletion, sum_pulse_runoff
  budget <- (st$total - start) - (flux[[2]] - flux[[3]] - flux[[4]])
  expect_gt(flux[[2]], 0)
  expect_gt(st$cover, 0)
  expect_lt(abs(budget), 1e-6 * flux[[2]])
})
