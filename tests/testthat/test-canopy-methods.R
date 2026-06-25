# Tests for the FF16 crown shading models (control$shading_model):
#   "deep-crown" - assimilation integrated over crown depth (default)
#   "crown-centre"   - assimilation a single evaluation at the crown centre; the
#                  light profile is built exactly as for deep-crown
#   "ppa"        - as crown-centre for assimilation, but the patch light profile is
#                  built as a stepped (layered) function of height
#
# All three share the same per-plant competition contribution (smooth Yokozawa
# Q). The model is resolved once in FF16_Strategy::prepare_strategy() into a
# function pointer (assimilation_fn) and, for the profile, configured on the
# environment by the Patch constructor -- so it costs no per-call string
# comparison on the hot path.
#
# Most SCM-level tests here run on a shortened patch horizon
# (max_patch_lifetime = 40 vs the FF16 default ~105) for speed: their assertions
# are relational (one model > another, or two models differ/agree at the same
# horizon), which the shorter patch preserves. The single exception is the
# "deep-crown reproduces the baseline SCM result" test, kept at the full default
# horizon so it still anchors the canonical full-lifetime FF16 number.
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

models <- c("deep-crown", "mean-light", "crown-centre", "ppa")

test_that("control defaults", {
  # "" means "use the strategy's own default" (FF16 -> deep-crown,
  # TF24 -> mean-light).
  expect_equal(Control()$shading_model, "")
  expect_equal(Control()$ppa_layer_optical_depth, 0.5)
  expect_equal(Control()$ppa_layer_smoothing, 0.3)
})

test_that("unknown shading model is rejected at strategy preparation", {
  s <- FF16_Strategy()
  s$control$shading_model <- "not-a-model"
  expect_error(FF16_Individual(s), "Unknown shading model: not-a-model")
})

test_that("every FF16 shading model prepares and computes at the individual level", {
  # All six models, including the box variants (whose competition differs but
  # which still compute a finite individual carbon balance).
  for (m in c(models, "flat-top-box", "flat-top-soft-box")) {
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
  for (m in c("mean-light", "crown-centre", "ppa")) {
    other <- sapply(zs, function(z) make_ind(m, h)$compute_competition(z))
    expect_equal(other, ref, tolerance = 1e-12)
  }
  # and it is a smooth, monotonically declining profile (not a step)
  expect_true(all(diff(ref) <= 1e-12))
  expect_equal(tail(ref, 1), 0)
})

test_that("flat-top-box casts a step competition profile, unlike crown-centre", {
  h <- 10
  eta_c <- local({ eta <- FF16_Strategy()$pars$eta; 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta) })
  zs <- seq(0, h, length.out = 41)

  smooth <- sapply(zs, function(z) make_ind("crown-centre", h)$compute_competition(z))
  box    <- sapply(zs, function(z) make_ind("flat-top-box", h)$compute_competition(z))

  # crown-centre is the smooth Yokozawa profile (many distinct interior values)
  expect_gt(length(unique(round(smooth[zs > 0 & zs < h], 6))), 2)
  # flat-top-box is a step: a single full value below the crown centre, 0 above
  full <- box[[1]]
  expect_equal(box[zs < h * eta_c], rep(full, sum(zs < h * eta_c)))
  expect_equal(box[zs > h * eta_c], rep(0, sum(zs > h * eta_c)))
})

test_that("flat-top-soft-box has a continuous competition profile and runs", {
  h <- 10
  zs <- seq(0, h, length.out = 200)
  box  <- sapply(zs, function(z) make_ind("flat-top-box", h)$compute_competition(z))
  soft <- sapply(zs, function(z) make_ind("flat-top-soft-box", h)$compute_competition(z))

  # The soft box removes the jump: its largest step between adjacent points is
  # far smaller than the hard box's (which drops the full value in one step).
  expect_lt(max(abs(diff(soft))), 0.25 * max(abs(diff(box))))
  expect_true(all(diff(soft) <= 1e-9))            # still monotone decreasing
  expect_equal(tail(soft, 1), 0)                  # zero leaf above the crown top

  # Being continuous, it builds a light environment and runs (unlike flat-top-box)
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 40
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  ctrl <- Control(); ctrl$shading_model <- "flat-top-soft-box"
  out <- run_scm(p1, Environment("FF16"), ctrl)
  expect_true(is.finite(out$offspring_production))

  # The mis-shaped competition propagates to the dynamics: a soft-box stand
  # differs from the correct crown-centre stand (same crown-centre assimilation, but
  # a different shade profile).
  ctrl_flat <- Control(); ctrl_flat$shading_model <- "crown-centre"
  out_flat <- run_scm(p1, Environment("FF16"), ctrl_flat)
  expect_false(isTRUE(all.equal(out$offspring_production,
                                out_flat$offspring_production, tolerance = 1e-2)))
})

test_that("flat-top-box cannot build a light environment (discontinuous competition)", {
  # The step competition makes the patch light profile discontinuous, so the
  # adaptive light-environment spline cannot represent it and the SCM fails.
  # This is the point of the model: the competition profile must be continuous.
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 40
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  ctrl <- Control(); ctrl$shading_model <- "flat-top-box"
  expect_error(run_scm(p1, Environment("FF16"), ctrl),
               "Interpolated function as refined as currently possible")
})

test_that("under uniform light, the integrate-based models all agree", {
  # With light constant in height, the leaf-area-weighted mean light equals the
  # light everywhere, and the crown leaf-density profile integrates to one. So
  # integrating photosynthesis over depth (deep-crown), integrating light then
  # evaluating once (mean-light), and a single evaluation at the crown centre
  # (crown-centre) all collapse to the same value.
  prod <- function(model, E) {
    ind <- make_ind(model)
    env <- Environment("FF16")
    env$set_fixed_environment(E, 100)
    ind$compute_rates(env)
    ind$aux("net_mass_production_dt")
  }
  for (E in c(1.0, 0.5, 0.2)) {
    ref <- prod("deep-crown", E)
    expect_equal(prod("mean-light", E), ref, tolerance = 1e-10)
    expect_equal(prod("crown-centre", E), ref, tolerance = 1e-10)
  }
})

test_that("mean-light assimilation >= deep-crown (Jensen, concave photosynthesis)", {
  # Photosynthesis saturates (is concave) in light, so evaluating it at the mean
  # light (mean-light) is >= the mean of the rate over the light distribution
  # (deep-crown). The two are equal only under uniform light.
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 40
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  run_op <- function(m) {
    ctrl <- Control(); ctrl$shading_model <- m
    run_scm(p1, Environment("FF16"), ctrl)$offspring_production
  }
  op_deep <- run_op("deep-crown")
  op_avg  <- run_op("mean-light")
  expect_true(is.finite(op_avg))
  expect_gt(op_avg, op_deep)
})

test_that("deep-crown reproduces the baseline SCM result", {
  # The default model must be the established FF16 behaviour. This is the one
  # SCM test kept at the full default horizon: it is the canonical anchor, so it
  # pins the established full-lifetime FF16 number rather than a shortened-horizon
  # value.
  p0 <- scm_base_parameters("FF16")
  env <- Environment("FF16")
  ctrl <- Control() # shading_model defaults to "deep-crown"
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  out <- run_scm(p1, env, ctrl)
  expect_equal(out$offspring_production, 16.88946, tolerance = 1e-4)
})

test_that("crown-centre runs through the SCM and changes the outcome", {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 40
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  # Compare crown-centre against deep-crown at the *same* (shortened) horizon, so
  # any difference is purely the model, not the horizon -- no dependence on the
  # full-lifetime baseline constant.
  deep <- run_scm(p1, Environment("FF16"), Control())$offspring_production
  ctrl <- Control(); ctrl$shading_model <- "crown-centre"
  out <- run_scm(p1, Environment("FF16"), ctrl)
  expect_true(is.finite(out$offspring_production))
  # crown-centre removes self-shading within the crown, so production differs
  expect_false(isTRUE(all.equal(out$offspring_production, deep, tolerance = 1e-3)))
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

  # crown-centre and deep-crown leave the profile smooth
  for (m in c("crown-centre", "deep-crown")) {
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
  p0$max_patch_lifetime <- 40
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  ctrl_deep <- Control(); ctrl_deep$shading_model <- "deep-crown"
  ctrl_ppa  <- Control(); ctrl_ppa$shading_model  <- "ppa"
  out_deep <- run_scm(p1, Environment("FF16"), ctrl_deep)
  out_ppa  <- run_scm(p1, Environment("FF16"), ctrl_ppa)
  expect_true(is.finite(out_ppa$offspring_production))
  expect_gt(out_ppa$offspring_production, out_deep$offspring_production)
})

test_that("a hard PPA step (no smoothing) defeats the adaptive solver", {
  # ppa_layer_smoothing = 0 is a genuine discontinuity, which the error-controlled
  # solver cannot integrate (it shrinks the step indefinitely). This is why the
  # default smoothing exists.
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 40
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  ctrl <- Control(); ctrl$shading_model <- "ppa"; ctrl$ppa_layer_smoothing <- 0
  expect_error(run_scm(p1, Environment("FF16"), ctrl))
})

test_that("adaptive and fixed-schedule PPA agree (well-behaved integration)", {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 40
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  ctrl <- Control(); ctrl$shading_model <- "ppa"
  adaptive <- run_scm(p1, Environment("FF16"), ctrl)$offspring_production
  pf <- p1; pf$ode_times <- seq(0, p0$max_patch_lifetime, length.out = 800)
  fixed <- run_scm(pf, Environment("FF16"), ctrl,
                   use_ode_times = TRUE)$offspring_production
  expect_equal(adaptive, fixed, tolerance = 1e-2)
})

# ---------------------------------------------------------------------------
# TF24 supports deep-crown, mean-light (its default) and crown-centre, but not
# PPA. The shading model controls how the (expensive) hydraulic leaf
# optimisation is aggregated over the crown: one evaluation at the mean light
# (mean-light) or crown-centre light (crown-centre), or one per crown-depth
# quadrature point with all leaf outputs integrated (deep-crown).
# ---------------------------------------------------------------------------

tf24_ind <- function(model, height = 5) {
  s <- TF24_Strategy()
  if (nzchar(model)) s$control$shading_model <- model
  ind <- TF24_Individual(s)
  ind$set_state("height", height)
  ind
}

tf24_prod <- function(model, E = 0.6) {
  ind <- tf24_ind(model)
  env <- Environment("TF24")
  env$set_fixed_environment(E, 50)
  ind$compute_rates(env)
  ind$aux("net_mass_production_dt")
}

test_that("TF24 defaults to mean-light and rejects PPA", {
  # The empty Control default maps to TF24's own default, mean-light, so
  # default behaviour is unchanged.
  expect_equal(tf24_prod(""), tf24_prod("mean-light"), tolerance = 1e-12)

  for (m in c("ppa", "flat-top-box", "flat-top-soft-box")) {
    s <- TF24_Strategy()
    s$control$shading_model <- m
    expect_error(TF24_Individual(s), "not supported for the TF24 strategy")
  }
})

test_that("TF24 shading models agree under uniform light", {
  # With light constant in height, the crown-centre light, the leaf-area-weighted
  # mean light, and the depth integral of the leaf optimisation all coincide.
  for (E in c(0.4, 0.7, 1.0)) {
    ref <- tf24_prod("mean-light", E)
    expect_equal(tf24_prod("crown-centre", E), ref, tolerance = 1e-8)
    expect_equal(tf24_prod("deep-crown", E), ref, tolerance = 1e-8)
  }
})
