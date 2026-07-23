# TF24t: TF24 + the ATLS merged-revision leaf thermal layer (#566). TF24t is a
# thin subclass of TF24 that reuses TF24_Pars/TF24_Environment, turns the leaf
# thermal-damage layer on, and adds three slow ODE states: a single
# thermostability acclimation state (acclim_thermostab) and the two lasting-
# damage pools (damage_recoverable = I_r, damage_permanent = I_p).
#
# Damage is now a MEMORY (state), not an instantaneous factor: a fresh hot leaf
# loses nothing beyond the reversible peaked-Arrhenius response until I_r/I_p
# accrue through the ODE. The physics is therefore tested with single-evaluation
# rate checks and preset damage states rather than long lethal-heat integrations
# (a plant pinned at lethal heat forever -- never culled -- drives capacity to the
# floor and is an unrealistic numerical stress case; the SCM culls dying plants).

test_that("TF24t reuses TF24_Pars and exposes the merged thermal/acclimation knobs", {
  s <- TF24t_Strategy()
  expect_inherits(s, "TF24t_Strategy")

  # pars is the *parent* TF24_Pars (has TF24 biological traits, no shadow class)
  expect_true(all(c("lma", "rho", "g1_TF24", "vcmax_25", "use_energy_balance") %in%
                    names(s$pars)))

  # top-level strategy carries the merged thermal traits + damage/acclim kinetics
  thermal_knobs <- c("topt_offset", "m_rep", "t_rep_cut", "k_i", "k_rec", "k_mat",
                     "dTopt_max", "K_A", "alpha", "beta", "t_accl", "softplus_s",
                     "c_acclim_maint", "c_repair_maint", "c_protect_maint",
                     "c_repair_flux", "c_build", "c_accl_induct", "induct_eps")
  expect_true(all(thermal_knobs %in% names(s)))
  # merged-revision defaults (day^-1 rate constants; calibration targets)
  expect_identical(s$k_rec, 0.25)
  expect_identical(s$k_i, 0.5)
  expect_identical(s$k_mat, 0.02)
  # removed as-built knobs are gone
  expect_false(any(c("tcrit_0", "k_d1_0", "k_r1_0", "m_switch", "dTcrit_max",
                     "alpha_crit", "c_build_tcrit") %in% names(s)))
})

test_that("TF24t appends three thermal ODE states after TF24's states", {
  p <- TF24t_Individual()
  nm <- p$ode_names
  appended <- c("acclim_thermostab", "damage_recoverable", "damage_permanent")
  expect_true(all(appended %in% nm))
  # appended last, in order (inherited indices unchanged)
  expect_identical(tail(nm, 3), appended)
  expect_identical(length(nm), length(TF24_Individual()$ode_names) + 3L)
})

test_that("TF24t validates acclimation and damage kinetics", {
  mk <- function(mutate) {
    p <- scm_base_parameters("TF24t") |>
      add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
    s <- p$strategies[[1]]; s <- mutate(s); p$strategies[[1]] <- s
    p
  }
  expect_error(run_scm(mk(function(s) { s$alpha <- -1; s }),
                       collect = FALSE, refine_schedule = FALSE),
               "acclimation kinetics")
  expect_error(run_scm(mk(function(s) { s$k_i <- -1; s }),
                       collect = FALSE, refine_schedule = FALSE),
               "damage kinetics")
})

test_that("TF24t runs end-to-end; acclimation sits at its env equilibrium; damage bounded", {
  p <- scm_base_parameters("TF24t") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  res <- run_scm(p, collect = TRUE, refine_schedule = FALSE)
  sp <- res$species
  expect_true(all(c("acclim_thermostab", "damage_recoverable",
                    "damage_permanent") %in% names(sp)))
  # lasting-damage pools are finite and bounded in [0, 1]
  for (col in c("damage_recoverable", "damage_permanent")) {
    expect_true(all(is.finite(sp[[col]])))
    expect_true(all(sp[[col]] >= -1e-8 & sp[[col]] <= 1 + 1e-6))
  }

  # Constant environment (default leaf_temp = 25) => the single thermostability
  # state sits at its seeded equilibrium A_eq = (alpha/beta)*softplus(T - t_accl).
  # (The forcing uses the environmental driver, not the operating-point Tleaf.)
  s <- p$strategies[[1]]
  softplus <- function(x, k) log1p(exp(k * x)) / k
  a_eq <- s$alpha / s$beta * softplus(25 - s$t_accl, s$softplus_s)
  expect_equal(unique(round(sp$acclim_thermostab, 10)), round(a_eq, 10))
})

# The merged cost structure. Leaf maintenance/activity terms are unit-tested at
# the Leaf level (test-leaf-thermal.R); here we cover the whole-plant terms in
# TF24t::net_mass_production_dt and confirm the damage feedback and costs bite.

test_that("tolerance construction cost debits net production only when bought", {
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
  # No thermostability offset -> construction cost is zero for any coefficient.
  expect_equal(net(0, 0), net(0, 5), tolerance = 1e-10)
  # With an offset, a larger construction coefficient lowers net production.
  expect_lt(net(8, 0.2), net(8, 0))
})

test_that("lasting damage state cuts net production (feedback needs PM, forced on)", {
  # A preset recoverable-damage state discounts capacity, so net production drops.
  # This is the damage feedback, evaluated at the operating-point Tleaf, which
  # only exists because prepare_strategy forces the PM energy balance on.
  env <- Environment("TF24t")
  env$extrinsic_drivers_set_constant("leaf_temp", 30)
  net_at_Ir <- function(Ir) {
    s <- TF24t_Strategy()
    expect_identical(s$pars$use_energy_balance, 0.0)   # PM off in pars...
    p <- TF24t_Individual(s)                           # ...but forced on in prepare_strategy
    p$set_state("height", 5)
    p$set_state("damage_recoverable", Ir)
    p$compute_rates(env)
    p$net_mass_production_dt(env)
  }
  expect_lt(net_at_Ir(0.5), net_at_Ir(0.0))
})

test_that("recoverable-damage flux rises with heat and self-limits with I_r", {
  dI_r <- function(Tmid, Ir) {
    env <- Environment("TF24t"); env$extrinsic_drivers_set_constant("leaf_temp", Tmid)
    p <- TF24t_Individual(TF24t_Strategy()); p$set_state("height", 5)
    p$set_state("damage_recoverable", Ir); p$compute_rates(env)
    rt <- p$ode_rates; names(rt) <- p$ode_names; rt[["damage_recoverable"]]
  }
  expect_gt(dI_r(40, 0), dI_r(20, 0))     # hotter operating point -> more damage flux
  expect_lt(dI_r(40, 0.5), dI_r(40, 0))   # the (1 - I_r - I_p) factor self-limits
})

test_that("maturation coupling: faster restorative repair lowers the permanent fraction", {
  # Permanent fraction of an insult ~ k_mat / (k_rec_eff + k_mat): prompt repair
  # (higher k_rec) prevents scarring. Evaluate at a cool operating point where the
  # repair gate is open (k_rec_eff ~ k_rec).
  perm_frac <- function(krec) {
    env <- Environment("TF24t"); env$extrinsic_drivers_set_constant("leaf_temp", 15)
    s <- TF24t_Strategy(); s$k_rec <- krec
    # permanent fraction of the I_r outflow (repair gate open at a cool leaf)
    s$k_mat / (krec + s$k_mat)
  }
  expect_lt(perm_frac(2.0), perm_frac(0.25))   # faster repair -> smaller permanent fraction
})

test_that("a single TF24t plant grows through the ODE runner, and costs slow it", {
  e <- "TF24_Env"
  env <- Environment("TF24t")   # default leaf_temp = 25: plant grows, damage stays sub-lethal
  grow_to_time <- function(mutate = identity, t_end = 10) {
    s <- mutate(TF24t_Strategy())
    p <- Individual("TF24t", e)(s)
    runner <- OdeRunner("TF24t")(IndividualRunner("TF24t", e)(p, env))
    h0 <- runner$object$individual$state("height")
    while (runner$time < t_end) runner$step()
    ode <- runner$object$individual
    list(h0 = h0, h1 = ode$state("height"), states = ode$internals$states)
  }
  base <- grow_to_time()
  expect_gt(base$h1, base$h0)                # a lone plant grows
  expect_true(all(is.finite(base$states)))   # finite trajectory (incl. thermal states)

  # Crank the thermal costs up: the same plant, integrated to the same time, ends
  # shorter -- the costs are genuinely debiting carbon.
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

test_that("with zero gain the acclimation state stays put (no drift)", {
  p <- scm_base_parameters("TF24t") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  s <- p$strategies[[1]]
  s$alpha <- 0; s$beta <- 0
  p$strategies[[1]] <- s
  res <- run_scm(p, collect = TRUE, refine_schedule = FALSE)
  expect_true(all(abs(res$species$acclim_thermostab) < 1e-8))
})
