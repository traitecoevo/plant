# ES20: FF16 with a non-structural carbohydrate store (Stefaniak et al.,
# egusphere-2026-1474). See inst/include/plant/models/es20_model.h.

es20_strategy <- function(t_s = 0.5, a_s = 0.1 * 365) {
  s <- ES20_Model()
  s$pars$t_s <- t_s
  s$pars$a_s <- a_s
  s
}

# One ES20 plant grown alone in full light.
es20_grow <- function(s, onset, times) {
  env <- ES20_Environment(onset)
  env$set_fixed_environment(1.0, 150.0)
  run(ES20_Individual(s, env), env = env, times = times)
}

test_that("ES20 declares its states and runs in FF16's environment", {
  s <- ES20_Model()
  expect_identical(environment_type("ES20"), "FF16_Env")
  expect_identical(model_id("ES20"), "ES20@v1")
  pl <- ES20_Individual(s)
  expect_identical(pl$ode_names[1:6],
                   c("height", "mortality", "fecundity", "area_heartwood",
                     "mass_heartwood", "mass_storage"))
  expect_length(pl$ode_names, 15L)
})

test_that("a new ES20 plant is seeded with its tissues and storage pool", {
  s <- es20_strategy()
  pl <- ES20_Individual(s)
  live <- pl$state("mass_leaf") + pl$state("mass_sapwood") +
    pl$state("mass_bark") + pl$state("mass_root")
  # The seedling's live mass is the seed mass, and storage is b_s1 of it (Eq. 8).
  expect_equal(live, s$pars$omega, tolerance = 1e-6)
  expect_equal(pl$state("mass_storage"), s$pars$b_s1 * live)
  expect_equal(pl$state("mass_leaf"), s$pars$lma * pl$state("area_leaf"))
  expect_equal(pl$aux("competition_effect"), pl$state("area_leaf"))
})

test_that("the stress factor follows the stress_onset driver", {
  s <- es20_strategy()
  stress_at <- function(env, time) {
    env$set_fixed_environment(1.0, 150.0)
    pl <- ES20_Individual(s, env)
    env$time <- time
    pl$compute_rates(env)
    pl$aux("stress")
  }
  # No driver: no stress, at any time of year.
  expect_identical(stress_at(ES20_Environment(), 0.9), 1.0)
  # Onset 0.75: unstressed mid-growing season, stressed 0.1 yr into the stress
  # period (the transitions are logistic, about a tenth of a year wide).
  env <- ES20_Environment(rep(0.75, 5))
  expect_gt(stress_at(env, 2.4), 1 - 1e-6)
  expect_lt(stress_at(env, 2.85), 1e-2)
  # A later onset that year leaves the plant unstressed at the same time.
  expect_gt(stress_at(ES20_Environment(c(0.75, 0.75, 0.95, 0.75)), 2.85), 1 - 1e-2)
})

test_that("stress halts photosynthesis and establishment and cuts respiration", {
  s <- es20_strategy()
  env <- ES20_Environment(rep(0.75, 5))
  env$set_fixed_environment(1.0, 150.0)
  pl <- ES20_Individual(s, env)
  env$time <- 0.4
  pl$compute_rates(env)
  growing <- c(npp = pl$aux("net_mass_production_dt"),
               resp = pl$aux("respiration_dt"))
  expect_gt(growing[["npp"]], 0)
  expect_gt(pl$establishment_probability(env), 0)
  env$time <- 0.875
  pl$compute_rates(env)
  # Assimilation is ~0, so net production is ~minus the (halved) respiration.
  expect_lt(pl$aux("net_mass_production_dt"), 0)
  # The stress factor is ~0.002 here, not 0, so these hold to about 1%.
  expect_equal(pl$aux("respiration_dt") / growing[["resp"]],
               s$pars$r_min, tolerance = 1e-2)
  expect_equal(pl$aux("net_mass_production_dt"), -pl$aux("respiration_dt"),
               tolerance = 2e-2)
  expect_identical(pl$establishment_probability(env), 0.0)
})

test_that("growth stops after the switch time and storage builds", {
  s <- es20_strategy(t_s = 0.25)
  res <- es20_grow(s, NULL, seq(0.02, 1, by = 0.02))
  h <- res$state[, "height"]
  st <- res$state[, "mass_storage"]
  t <- res$time
  at <- function(x) which.min(abs(t - x))
  # Height grows before t_s and barely moves between t_s + 0.1 and the year end.
  expect_gt(h[at(0.2)] - h[at(0.02)], 0)
  late <- t > 0.35 & t < 0.95
  expect_lt(diff(range(h[late])) / h[at(0.2)], 1e-3)
  expect_gt(st[at(0.94)], st[at(0.36)])
})

test_that("a Safe strategy ends a stressed year with more storage than a Risky one", {
  onset <- rep(0.75, 5)
  times <- seq(0.05, 2, by = 0.05)
  safe <- es20_grow(es20_strategy(t_s = 0.25), onset, times)
  risky <- es20_grow(es20_strategy(t_s = 0.5), onset, times)
  at <- which.min(abs(times - 1.7))
  portion <- function(r) {
    x <- r$state[at, ]
    x[["mass_storage"]] / (x[["mass_leaf"]] + x[["mass_sapwood"]] +
                             x[["mass_bark"]] + x[["mass_root"]])
  }
  expect_gt(portion(safe), portion(risky))
  expect_gt(risky$state[at, "height"], safe$state[at, "height"])
})

test_that("ES20 runs a patch by both solvers", {
  p <- base_parameters("ES20")
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(c(0.25, 0.5), "t_s"), birth_rate = c(0.5, 0.5))
  env <- ES20_Environment(rep(0.75, 10))
  det <- run(p, env = env)
  expect_s3_class(det, "plant_run")
  set.seed(1)
  p$patch_area <- 100
  sto <- run(p, env = env, method = "stochastic")
  expect_s3_class(sto, "plant_run")
})
