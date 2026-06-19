# Tests for the FF16 crown shading models (control$shading_model):
#   "deep-crown" - assimilation integrated over crown depth (default)
#   "flat-top"   - assimilation a single evaluation at the crown centre; the
#                  light profile is built exactly as for deep-crown
#   "ppa"        - as flat-top for assimilation, but the patch light profile is
#                  built as a stepped (layered) function of height
#
# All three share the same per-plant competition contribution (smooth Yokozawa
# Q). The model is resolved once in FF16_Strategy::prepare_strategy() into a
# function pointer (assimilation_fn) and, for the profile, configured on the
# environment by the Patch constructor -- so it costs no per-call string
# comparison on the hot path.
context("Canopy shading methods")

# A prepared FF16 individual under a given shading model. Constructing the
# Individual triggers make_strategy_ptr() -> prepare_strategy(), which binds the
# assimilation function pointer. NOTE: a stand-alone Individual is not inside a
# Patch, so its environment is never put into the stepped (PPA) mode; PPA-vs-
# others differences in the light profile only appear through run_scm() below.
make_ind <- function(model, height = 10) {
  s <- FF16_Strategy()
  s$control$shading_model <- model
  ind <- FF16_Individual(s)
  ind$set_state("height", height)
  ind
}

models <- c("deep-crown", "flat-top", "ppa")

test_that("control defaults", {
  expect_equal(Control()$shading_model, "deep-crown")
  expect_equal(Control()$ppa_layer_optical_depth, 0.5)
  expect_equal(Control()$ppa_layer_smoothing, 0.3)
})

test_that("unknown shading model is rejected at strategy preparation", {
  s <- FF16_Strategy()
  s$control$shading_model <- "not-a-model"
  expect_error(FF16_Individual(s), "Unknown shading model: not-a-model")
})

test_that("all three models prepare and compute without error", {
  for (m in models) {
    ind <- make_ind(m)
    env <- Environment("FF16")
    env$set_fixed_environment(0.5, 100)
    expect_silent(ind$compute_rates(env))
    expect_true(is.finite(ind$aux("net_mass_production_dt")))
  }
})

test_that("per-plant competition is identical across models (all use smooth Q)", {
  h <- 10
  zs <- seq(0, h, length.out = 21)
  ref <- sapply(zs, function(z) make_ind("deep-crown", h)$compute_competition(z))
  for (m in c("flat-top", "ppa")) {
    other <- sapply(zs, function(z) make_ind(m, h)$compute_competition(z))
    expect_equal(other, ref, tolerance = 1e-12)
  }
  # and it is a smooth, monotonically declining profile (not a step)
  expect_true(all(diff(ref) <= 1e-12))
  expect_equal(tail(ref, 1), 0)
})

test_that("under uniform light, deep-crown and flat-top assimilate identically", {
  # With light constant in height, integrating photosynthesis * leaf density over
  # crown depth (deep-crown) reduces exactly to a single evaluation at the crown
  # centre (flat-top), because the crown leaf-density profile integrates to one.
  for (E in c(1.0, 0.5, 0.2)) {
    deep <- make_ind("deep-crown")
    flat <- make_ind("flat-top")
    for (ind in list(deep, flat)) {
      env <- Environment("FF16")
      env$set_fixed_environment(E, 100)
      ind$compute_rates(env)
    }
    expect_equal(deep$aux("net_mass_production_dt"),
                 flat$aux("net_mass_production_dt"),
                 tolerance = 1e-10)
  }
})

test_that("deep-crown reproduces the baseline SCM result", {
  # The default model must be the established FF16 behaviour.
  p0 <- scm_base_parameters("FF16")
  env <- Environment("FF16")
  ctrl <- Control() # shading_model defaults to "deep-crown"
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16_hyperpar,
                          birth_rate_list = list(20))
  out <- run_scm(p1, env, ctrl)
  expect_equal(out$offspring_production, 16.88946, tolerance = 1e-4)
})

test_that("flat-top runs through the SCM and changes the outcome", {
  p0 <- scm_base_parameters("FF16")
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16_hyperpar,
                          birth_rate_list = list(20))
  ctrl <- Control(); ctrl$shading_model <- "flat-top"
  out <- run_scm(p1, Environment("FF16"), ctrl)
  expect_true(is.finite(out$offspring_production))
  # flat-top removes self-shading within the crown, so production differs
  expect_false(isTRUE(all.equal(out$offspring_production, 16.88946,
                                tolerance = 1e-3)))
})

test_that("PPA discretises the light profile into optical-depth layers", {
  # Directly exercise the stepped-profile transform (the PPA-specific machinery)
  # via set_shading_model(). With a hard step (smoothing -> 0) and a smooth light
  # value E, PPA returns exp(-d * floor(-log(E) / d)).
  hard <- 0.0 # layer_smoothing -> hard floor, for exact layer values
  env <- Environment("FF16")
  env$set_fixed_environment(0.3, 100)
  expect_equal(env$get_environment_at_height(10), 0.3) # smooth by default

  d <- 0.5
  env$set_shading_model("ppa", d, hard)
  # -log(0.3) = 1.20397; floor(1.20397 / 0.5) = 2; exp(-0.5 * 2) = exp(-1)
  expect_equal(env$get_environment_at_height(10), exp(-1.0), tolerance = 1e-9)

  # full light is unchanged (top of canopy: zero optical depth, zero layers)
  env$set_fixed_environment(1.0, 100)
  env$set_shading_model("ppa", d, hard)
  expect_equal(env$get_environment_at_height(10), 1.0)

  # coarser layers -> different step value: floor(1.20397 / 1.0) = 1 -> exp(-1)
  env$set_fixed_environment(0.3, 100)
  env$set_shading_model("ppa", 1.0, hard)
  expect_equal(env$get_environment_at_height(10), exp(-1.0), tolerance = 1e-9)

  # flat-top and deep-crown leave the profile smooth
  for (m in c("flat-top", "deep-crown")) {
    env$set_fixed_environment(0.3, 100)
    env$set_shading_model(m, d, hard)
    expect_equal(env$get_environment_at_height(10), 0.3)
  }
})

test_that("PPA layer smoothing keeps the profile monotone and bounded", {
  # The smoothed staircase must be monotone in height and stay within (0, 1].
  env <- Environment("FF16")
  env$set_fixed_environment(0.3, 100)
  env$set_shading_model("ppa", 0.5, 0.3) # default-style smoothing
  # set_fixed_environment is uniform, so vary E by sweeping the *value*
  Es <- seq(0.01, 1.0, length.out = 50)
  stepped <- sapply(Es, function(e) {
    env$set_fixed_environment(e, 100)
    env$set_shading_model("ppa", 0.5, 0.3)
    env$get_environment_at_height(10)
  })
  expect_true(all(stepped > 0 & stepped <= 1 + 1e-12))
  expect_true(all(diff(stepped) >= -1e-9)) # monotone increasing in E
})

test_that("PPA runs through the SCM (smoothed) and changes the outcome", {
  # With the default C1 smoothing (ppa_layer_smoothing = 0.3) the stepped profile
  # is differentiable, so the normal adaptive solver integrates PPA without a
  # fixed schedule. The layered light reduces self-shading, so production differs
  # markedly from deep-crown.
  p0 <- scm_base_parameters("FF16")
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16_hyperpar,
                          birth_rate_list = list(20))
  ctrl_deep <- Control(); ctrl_deep$shading_model <- "deep-crown"
  ctrl_ppa  <- Control(); ctrl_ppa$shading_model  <- "ppa"
  out_deep <- run_scm(p1, Environment("FF16"), ctrl_deep)
  out_ppa  <- run_scm(p1, Environment("FF16"), ctrl_ppa)
  expect_true(is.finite(out_ppa$offspring_production))
  expect_gt(out_ppa$offspring_production, out_deep$offspring_production)
})

test_that("adaptive and fixed-schedule PPA agree (well-behaved integration)", {
  p0 <- scm_base_parameters("FF16")
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16_hyperpar,
                          birth_rate_list = list(20))
  ctrl <- Control(); ctrl$shading_model <- "ppa"
  adaptive <- run_scm(p1, Environment("FF16"), ctrl)$offspring_production
  pf <- p1; pf$ode_times <- seq(0, p0$max_patch_lifetime, length.out = 2000)
  fixed <- run_scm(pf, Environment("FF16"), ctrl,
                   use_ode_times = TRUE)$offspring_production
  expect_equal(adaptive, fixed, tolerance = 1e-2)
})
