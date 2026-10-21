# TF24t: TF24 with phylloptim's lasting leaf thermal-damage and acclimation
# layer switched on (#566). TF24t reuses TF24_Pars and TF24_Environment, forces
# the leaf's energy balance on, and appends four ODE states: the thermostability
# (acclimation) state, the recoverable and permanent damaged fractions, and a
# running mean of the operating-point leaf temperature that feeds back into
# nothing.
#
# Damage is a STATE, not an instantaneous factor: a fresh hot leaf loses nothing
# beyond the reversible peaked Arrhenius until I_r and I_p accrue through the
# ODE. So the physics is tested on single rate evaluations at preset states, and
# the SCM runs are short.
#
# The leaf-level layer itself (the phi_d curve, the survival discount, the repair
# gate, the costs in R_d) is tested in phylloptim (test-thermal-damage.R).

tf24t_days <- 365.25

# phylloptim's helpers, in R, for the closed forms the rates are checked against.
tf24t_softplus <- function(x, s) log1p(exp(pmin(pmax(s * x, -30), 30))) / s
tf24t_logistic <- function(x) 1 / (1 + exp(-pmin(pmax(x, -30), 30)))

# The deactivated fraction at leaf temperature `Tleaf` (deg C) with the onset
# moved up by `delta` deg C: K/(1+K) of the Jmax curve's two-state equilibrium.
# H_d and d_S are phylloptim's Jmax defaults and the gas constant is its exact
# one; a change to either there moves this closed form.
tf24t_phi_d <- function(Tleaf, delta, H_d = 200000, d_S = 650) {
  R <- 8.314462618153240
  onset <- min(max(H_d / d_S + delta, 250), 350)
  T_K <- Tleaf + 273.15
  arg <- (H_d / onset * T_K - H_d) / (R * T_K)
  K <- exp(min(max(arg, -30), 30))
  K / (1 + K)
}

# Every thermal input at a setting that can neither damage nor cost: no
# acclimation gain (A stays at its zero seed), no leak into damage, no offset,
# and every cost coefficient zero.
tf24t_inert <- function(s = TF24t_Strategy()) {
  s$alpha <- 0
  s$k_i <- 0
  s$topt_offset <- 0
  for (k in c("c_acclim_maint", "c_repair_maint", "c_protect_maint",
              "c_repair_flux", "c_build", "c_accl_induct")) {
    s[[k]] <- 0
  }
  s
}

tf24t_rates <- function(ind) stats::setNames(ind$ode_rates, ind$ode_names)

test_that("Defaults", {
  s <- TF24t_Strategy()
  expect_inherits(s, "TF24t_Strategy")

  # pars is TF24's own, with no TF24t-specific shadow class.
  expect_identical(unclass(s$pars), unclass(TF24_Strategy()$pars))
  expect_identical(s$pars$use_energy_balance, 0.0)

  thermal <- list(topt_offset = 0, m_rep = 0.4, t_rep_cut = 45, k_i = 0.5,
                  k_rec = 0.25, k_mat = 0.02, dTopt_max = 6, K_A = 1,
                  alpha = 0.02, beta = 0.01, t_accl = 30, softplus_s = 1,
                  k_mean_leaf_temp = 0.0333,
                  c_acclim_maint = 0.02, c_repair_maint = 1e-5,
                  c_protect_maint = 1e-5, c_repair_flux = 1e-4,
                  c_build = 0.02, c_accl_induct = 0.01, induct_eps = 1e-4)
  expected_top <- c("pars", "control", "collect_all_auxiliary",
                    "birth_rate_x", "birth_rate_y", "is_variable_birth_rate",
                    names(thermal))
  expect_identical(sort(names(s)), sort(expected_top))
  expect_identical(unclass(s)[names(thermal)], thermal)
  expect_identical(s$control, Control())
})

test_that("TF24t appends four thermal ODE states after TF24's", {
  p <- TF24t_Individual()
  appended <- c("acclim_thermostab", "damage_recoverable", "damage_permanent",
                "mean_leaf_temp")
  expect_identical(p$ode_names, c(TF24_Individual()$ode_names, appended))
  # No aux of its own: the running mean relaxes toward TF24's `Tleaf`.
  expect_identical(p$aux_names, TF24_Individual()$aux_names)
})

test_that("the energy balance is forced on whatever pars says", {
  env <- Environment("TF24t")
  env$extrinsic_drivers_set_constant("leaf_temp", 30)
  s <- TF24t_Strategy()
  expect_identical(s$pars$use_energy_balance, 0.0)
  ind <- TF24t_Individual(s)
  ind$set_state("height", 5)
  ind$compute_rates(env)
  # Off the energy balance Tleaf IS the driver; a bright leaf solving its own
  # temperature sits elsewhere.
  expect_true(is.finite(ind$aux("Tleaf")))
  expect_false(isTRUE(all.equal(ind$aux("Tleaf"), 30)))

  ref <- TF24_Individual(TF24_Strategy())
  ref$set_state("height", 5)
  ref$compute_rates(env)
  expect_identical(ref$aux("Tleaf"), 30)
})

test_that("the thermal rates follow their equations at the operating point", {
  # Every trait the leaf reads is moved off its default, so a copy into the leaf
  # that prepare_strategy dropped shows up here as the leaf's default in phi_d or
  # k_rec_eff. phi_d and k_rec_eff are read at the operating-point leaf
  # temperature, which is TF24's `Tleaf` aux.
  air <- 33
  env <- Environment("TF24t")
  env$extrinsic_drivers_set_constant("leaf_temp", air)
  s <- TF24t_Strategy()
  s$topt_offset <- 1.5; s$dTopt_max <- 4; s$K_A <- 0.5
  s$k_rec <- 0.4; s$m_rep <- 0.3; s$t_rep_cut <- 40
  s$k_i <- 0.7; s$k_mat <- 0.05; s$alpha <- 0.03; s$beta <- 0.02
  s$k_mean_leaf_temp <- 0.05
  A <- 2; Ir <- 0.2; Ip <- 0.1; Tbar <- 28
  ind <- TF24t_Individual(s)
  ind$set_state("height", 5)
  ind$set_state("acclim_thermostab", A)
  ind$set_state("damage_recoverable", Ir)
  ind$set_state("damage_permanent", Ip)
  ind$set_state("mean_leaf_temp", Tbar)
  ind$compute_rates(env)
  rt <- tf24t_rates(ind)
  Tleaf <- ind$aux("Tleaf")
  expect_true(is.finite(Tleaf))

  # Acclimation is forced by the AIR temperature, the driver.
  expect_equal(rt[["acclim_thermostab"]],
               tf24t_days * (s$alpha * tf24t_softplus(air - s$t_accl, s$softplus_s) -
                             s$beta * A),
               tolerance = 1e-12)

  delta <- s$topt_offset + s$dTopt_max * A / (A + s$K_A)
  phi_d <- tf24t_phi_d(Tleaf, delta)
  k_rec_eff <- s$k_rec * tf24t_logistic(-s$m_rep * (Tleaf - s$t_rep_cut))
  expect_gt(phi_d, 0)
  expect_equal(rt[["damage_recoverable"]],
               tf24t_days * (s$k_i * phi_d * (1 - Ir - Ip) -
                             (k_rec_eff + s$k_mat) * Ir),
               tolerance = 1e-10)

  # Permanent damage is diluted by leaf turnover and canopy expansion, both
  # already per year; the expansion is recovered from the height rate through
  # TF24's allometry h = a_l1 * L^a_l2.
  pars <- s$pars
  L <- ind$aux("competition_effect")
  dL_dt <- rt[["height"]] / (pars$a_l1 * pars$a_l2 * L^(pars$a_l2 - 1))
  dilution <- pars$k_l + max(0, dL_dt) / L
  expect_equal(rt[["damage_permanent"]],
               tf24t_days * s$k_mat * Ir - dilution * Ip, tolerance = 1e-10)

  expect_equal(rt[["mean_leaf_temp"]],
               tf24t_days * s$k_mean_leaf_temp * (Tleaf - Tbar),
               tolerance = 1e-12)
})

test_that("each leaf cost coefficient reaches the leaf", {
  # The four maintenance costs are added to the leaf's R_d and so reach net
  # production only through the leaf. Each multiplicand is made positive -- A > 0,
  # I_r > 0, and k_i below the leaf's reference 0.5 so some protection is bought.
  env <- Environment("TF24t")
  env$extrinsic_drivers_set_constant("leaf_temp", 30)
  net <- function(mutate = identity) {
    s <- TF24t_Strategy()
    s$k_i <- 0.2
    s <- mutate(s)
    ind <- TF24t_Individual(s)
    ind$set_state("height", 5)
    ind$set_state("acclim_thermostab", 2)
    ind$set_state("damage_recoverable", 0.2)
    ind$compute_rates(env)
    ind$aux("net_mass_production_dt")
  }
  base <- net()
  expect_true(is.finite(base))
  for (k in c("c_acclim_maint", "c_repair_maint", "c_protect_maint",
              "c_repair_flux")) {
    expect_lt(net(function(s) { s[[k]] <- 10 * s[[k]] + 0.1; s }), base,
              label = k)
  }
  # k_i reaches the leaf too, through the protection it buys.
  expect_lt(net(function(s) { s$k_i <- 0; s$c_protect_maint <- 0.1; s }),
            net(function(s) { s$k_i <- 0.4; s$c_protect_maint <- 0.1; s }))
})

test_that("at non-damaging settings TF24t is TF24 with the energy balance on", {
  # What can be asserted, and why.
  #
  # With A held at zero, no leak into damage, no offset and every cost zero, the
  # layer multiplies both capacities by exactly 1 and adds exact zeros to R_d,
  # and the onset shift is zero. The one operation that is not an identity by
  # construction is the shift's round trip through the optimum, d_S -> T_opt ->
  # d_S, which is exact at phylloptim's default temperature parameters. So a rate
  # evaluation is bit-identical to TF24 run with pars$use_energy_balance = 1,
  # which is the energy balance TF24t forces.
  #
  # A free SCM run is NOT bit-identical, and should not be: TF24t's four extra
  # states sit in the adaptive step-size controller's error norm, so the two runs
  # take different steps (measured 260 against 266, offspring 7.7e-4 apart). What
  # does hold exactly is TF24 pinned to TF24t's own step program, which is the
  # assertion below. A change that leaks the layer into the shared states at
  # non-damaging settings breaks it; a step-controller difference cannot.
  env <- Environment("TF24t")
  for (air in c(25, 35)) {
    env$extrinsic_drivers_set_constant("leaf_temp", air)
    a <- TF24t_Individual(tf24t_inert())
    a$set_state("height", 5)
    a$compute_rates(env)
    s_ref <- TF24_Strategy()
    s_ref$pars$use_energy_balance <- 1
    expect_identical(s_ref$pars$use_energy_balance, 1)
    b <- TF24_Individual(s_ref)
    b$set_state("height", 5)
    b$compute_rates(env)
    shared <- b$ode_names
    expect_identical(tf24t_rates(a)[shared], tf24t_rates(b)[shared])
    expect_identical(a$internals$auxs, b$internals$auxs)
    expect_identical(unname(tf24t_rates(a)[c("acclim_thermostab",
                                             "damage_recoverable",
                                             "damage_permanent")]),
                     c(0, 0, 0))
  }

  build <- function(type, mutate) {
    p <- scm_base_parameters(type)
    p$max_patch_lifetime <- 5
    p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
    p$strategies[[1]] <- mutate(p$strategies[[1]])
    p
  }
  pm <- function(s) { s$pars$use_energy_balance <- 1; s }
  thermal <- run_scm(build("TF24t", tf24t_inert))
  p_ref <- build("TF24", pm)
  expect_identical(p_ref$strategies[[1]]$pars$use_energy_balance, 1)
  p_ref$ode_times <- thermal$ode_times
  p_ref$ode_step_sizes <- thermal$ode_step_sizes
  pinned <- run_scm(p_ref)

  expect_identical(thermal$offspring_production, pinned$offspring_production)
  expect_identical(thermal$net_reproduction_ratios,
                   pinned$net_reproduction_ratios)
  expect_identical(thermal$patch$environment$get_soil_water_state(),
                   pinned$patch$environment$get_soil_water_state())
  nodes_t <- thermal$patch$species[[1]]$nodes
  nodes_p <- pinned$patch$species[[1]]$nodes
  expect_identical(length(nodes_t), length(nodes_p))
  n_shared <- length(TF24_Individual()$ode_names)
  for (i in seq_along(nodes_t)) {
    expect_identical(nodes_t[[i]]$individual$ode_state[seq_len(n_shared)],
                     nodes_p[[i]]$individual$ode_state)
    expect_identical(nodes_t[[i]]$log_density, nodes_p[[i]]$log_density)
  }

  # Paired with the default TF24t, which the same comparison must tell apart.
  default <- run_scm(build("TF24t", identity))
  free_ref <- run_scm(build("TF24", pm))
  expect_gt(abs(log(default$offspring_production /
                    free_ref$offspring_production)), 1)
})

test_that("mean_leaf_temp is a running mean of the operating leaf temperature", {
  env <- Environment("TF24t")
  env$extrinsic_drivers_set_constant("leaf_temp", 30)

  # Seeded at the birth AIR temperature: the operating point is not solved until
  # the first compute_rates, which runs after set_initial_states.
  p0 <- TF24t_Individual(TF24t_Strategy())
  p0$set_initial_states(env)
  expect_identical(p0$state("mean_leaf_temp"), env$get_leaf_temp())

  rate_at <- function(mean0, k = 0.0333) {
    s <- TF24t_Strategy(); s$k_mean_leaf_temp <- k
    p <- TF24t_Individual(s); p$set_state("height", 5)
    p$set_state("mean_leaf_temp", mean0)
    p$compute_rates(env)
    list(rate = tf24t_rates(p)[["mean_leaf_temp"]], T_op = p$aux("Tleaf"))
  }
  hot <- rate_at(10); cold <- rate_at(60)
  expect_gt(hot$T_op, 30)   # a bright leaf heats above the 30 deg C air
  expect_gt(hot$rate, 0)    # below the operating temperature -> warming
  expect_lt(cold$rate, 0)   # above it -> cooling
  # The operating point does not read the running mean, so the rate is zero at
  # mean == T_op ...
  expect_equal(rate_at(hot$T_op)$rate, 0, tolerance = 1e-6)
  # ... and k = 0 freezes the state.
  expect_identical(rate_at(10, k = 0)$rate, 0)
})

test_that("TF24t refuses misconfigured kinetics", {
  mk <- function(mutate) {
    p <- scm_base_parameters("TF24t") |>
      add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
    p$strategies[[1]] <- mutate(p$strategies[[1]])
    p
  }
  refuses <- function(mutate, why) {
    expect_error(run_scm(mk(mutate), collect = FALSE, refine_schedule = FALSE),
                 why)
  }
  refuses(function(s) { s$alpha <- -1; s }, "acclimation kinetics")
  refuses(function(s) { s$beta <- -1; s }, "acclimation kinetics")
  refuses(function(s) { s$softplus_s <- 0; s }, "acclimation kinetics")
  refuses(function(s) { s$k_mean_leaf_temp <- -1; s }, "acclimation kinetics")
  refuses(function(s) { s$alpha <- NaN; s }, "acclimation kinetics")
  refuses(function(s) { s$k_i <- -1; s }, "damage kinetics")
  refuses(function(s) { s$k_rec <- -1; s }, "damage kinetics")
  refuses(function(s) { s$k_mat <- Inf; s }, "damage kinetics")
})

test_that("TF24t refuses deep-crown shading and accepts the other two", {
  # phi_d and k_rec_eff are read after the crown loop; deep-crown would leave
  # them at the last quadrature node.
  s <- TF24t_Strategy()
  s$control$shading_model <- "deep-crown"
  expect_error(TF24t_Individual(s), "deep-crown.*not supported for the TF24t")

  p <- scm_base_parameters("TF24t") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  ctrl <- control()
  ctrl$shading_model <- "deep-crown"
  expect_error(run_scm(p, ctrl = ctrl), "not supported for the TF24t")

  for (model in c("mean-light", "crown-centre")) {
    s <- TF24t_Strategy()
    s$control$shading_model <- model
    expect_inherits(TF24t_Individual(s), "Individual")
  }
})

test_that("TF24t has no stand gradient", {
  p <- scm_base_parameters("TF24t")
  p$max_patch_lifetime <- 1
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
  scm <- run_scm(p)
  expect_error(stand_gradient(scm), "thermal-damage layer")
})

test_that("TF24t runs end-to-end; acclimation sits at its equilibrium; damage is bounded", {
  p <- scm_base_parameters("TF24t")
  p$max_patch_lifetime <- 10
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
  res <- run_scm(p, collect = TRUE, refine_schedule = FALSE)
  sp <- res$species
  expect_true(all(c("acclim_thermostab", "damage_recoverable",
                    "damage_permanent", "mean_leaf_temp") %in% names(sp)))
  for (col in c("damage_recoverable", "damage_permanent")) {
    x <- sp[[col]][!is.na(sp[[col]])]
    expect_true(all(is.finite(x)), info = col)
    expect_true(all(x >= -1e-8 & x <= 1 + 1e-6), info = col)
  }
  # Damage does accrue at the defaults, so the bound above is not vacuous.
  expect_gt(max(sp$damage_recoverable, na.rm = TRUE), 0.01)

  # A constant environment (leaf_temp driver 25) keeps A at its seed,
  # A_eq = (alpha/beta) * softplus(T_air - t_accl).
  s <- p$strategies[[1]]
  a_eq <- s$alpha / s$beta * tf24t_softplus(25 - s$t_accl, s$softplus_s)
  a <- sp$acclim_thermostab[!is.na(sp$acclim_thermostab)]
  expect_equal(unique(round(a, 10)), round(a_eq, 10))
})

test_that("with zero gain the acclimation state stays at zero", {
  p <- scm_base_parameters("TF24t")
  p$max_patch_lifetime <- 10
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
  s <- p$strategies[[1]]
  s$alpha <- 0; s$beta <- 0
  p$strategies[[1]] <- s
  res <- run_scm(p, collect = TRUE, refine_schedule = FALSE)
  a <- res$species$acclim_thermostab[!is.na(res$species$acclim_thermostab)]
  expect_true(all(abs(a) < 1e-8))
})

# The whole-plant costs in net_mass_production_dt, and the damage feedback.

test_that("the construction cost debits net production only when bought", {
  env <- Environment("TF24t")
  net <- function(offset, cbuild) {
    s <- TF24t_Strategy()
    s$topt_offset <- offset
    s$c_build <- cbuild
    p <- TF24t_Individual(s)
    p$set_state("height", 5)
    p$compute_rates(env)
    p$net_mass_production_dt(env)
  }
  # No offset: the construction cost is zero whatever the coefficient.
  expect_identical(net(0, 0), net(0, 5))
  # With one, a larger coefficient lowers net production.
  expect_lt(net(8, 0.2), net(8, 0))
})

test_that("lasting damage cuts net production", {
  env <- Environment("TF24t")
  env$extrinsic_drivers_set_constant("leaf_temp", 30)
  net_at_Ir <- function(Ir) {
    p <- TF24t_Individual(TF24t_Strategy())
    p$set_state("height", 5)
    p$set_state("damage_recoverable", Ir)
    p$compute_rates(env)
    p$net_mass_production_dt(env)
  }
  expect_lt(net_at_Ir(0.5), net_at_Ir(0.0))
})

test_that("the recoverable-damage flux rises with heat and self-limits with I_r", {
  dI_r <- function(air, Ir) {
    env <- Environment("TF24t")
    env$extrinsic_drivers_set_constant("leaf_temp", air)
    p <- TF24t_Individual(TF24t_Strategy())
    p$set_state("height", 5)
    p$set_state("damage_recoverable", Ir)
    p$compute_rates(env)
    tf24t_rates(p)[["damage_recoverable"]]
  }
  expect_gt(dI_r(40, 0), dI_r(20, 0))     # hotter operating point, more flux
  expect_lt(dI_r(40, 0.5), dI_r(40, 0))   # the (1 - I_r - I_p) factor
})

test_that("faster repair lowers the share of damage that becomes permanent", {
  # With no inflow (k_i = 0) and no permanent damage yet, I_r drains at
  # (k_rec_eff + k_mat) and I_p gains k_mat of it, so the permanent share of the
  # outflow is k_mat / (k_rec_eff + k_mat). Read off the rates rather than
  # computed from the parameters, at a cool leaf where the repair gate is open.
  permanent_share <- function(k_rec) {
    env <- Environment("TF24t")
    env$extrinsic_drivers_set_constant("leaf_temp", 15)
    s <- TF24t_Strategy(); s$k_i <- 0; s$k_rec <- k_rec
    p <- TF24t_Individual(s)
    p$set_state("height", 5)
    p$set_state("damage_recoverable", 0.3)
    p$compute_rates(env)
    rt <- tf24t_rates(p)
    rt[["damage_permanent"]] / -rt[["damage_recoverable"]]
  }
  fast <- permanent_share(2.0)
  slow <- permanent_share(0.25)
  expect_lt(fast, slow)
  expect_gt(fast, 0)
})

test_that("a single TF24t plant grows through the ODE runner, and costs slow it", {
  # At 15 deg C air, not the default 25. A lone plant in full light runs its
  # leaves 6-9 deg C above the air, so at 25 phi_d sits near its midpoint, damage
  # ratchets into I_p within a year at the default k_mat, and the dying plant's
  # empty storage pool holds the stepper to tiny steps.
  #
  # The runner does not seed the initial states, so they are seeded here: an
  # unseeded storage pool starts empty, where a plant losing carbon is refused
  # every step that dips below it.
  e <- "TF24_Env"
  env <- Environment("TF24t")
  env$extrinsic_drivers_set_constant("leaf_temp", 15)
  grow_to_time <- function(mutate = identity, t_end = 10) {
    s <- mutate(TF24t_Strategy())
    p <- Individual("TF24t", e)(s)
    p$set_initial_states(env)
    runner <- OdeRunner("TF24t")(IndividualRunner("TF24t", e)(p, env))
    h0 <- runner$object$individual$state("height")
    while (runner$time < t_end) runner$step()
    ode <- runner$object$individual
    list(h0 = h0, h1 = ode$state("height"), states = ode$internals$states)
  }
  base <- grow_to_time()
  expect_gt(base$h1, base$h0)
  expect_true(all(is.finite(base$states)))

  costly <- grow_to_time(function(s) {
    s$topt_offset <- 8
    s$c_build <- 2
    s$c_acclim_maint <- 5; s$c_repair_maint <- 1e-2
    s$c_protect_maint <- 1e-2; s$c_repair_flux <- 1e-2
    s$c_accl_induct <- 5
    s
  })
  expect_true(all(is.finite(costly$states)))
  expect_lt(costly$h1, base$h1)
})
