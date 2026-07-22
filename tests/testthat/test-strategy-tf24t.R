# TF24t: TF24 + leaf thermal damage/repair/acclimation (ATLS, #566).
# TF24t is a thin subclass of TF24 that reuses TF24_Pars/TF24_Environment, turns
# the leaf thermal-damage layer on, and adds two acclimation ODE states
# (acclim_topt, acclim_tcrit).

test_that("TF24t reuses TF24_Pars and exposes the thermal/acclimation knobs", {
  s <- TF24t_Strategy()
  expect_inherits(s, "TF24t_Strategy")

  # pars is the *parent* TF24_Pars (has TF24 biological traits, no shadow class)
  expect_true(all(c("lma", "rho", "g1_TF24", "vcmax_25", "use_energy_balance") %in%
                    names(s$pars)))

  # top-level strategy carries the thermal-damage traits + acclimation kinetics
  thermal_knobs <- c("topt_offset", "tcrit_0", "k_d1_0", "k_r1_0", "m_switch",
                     "m_rep", "t_rep_cut", "dTcrit_max", "dTopt_max", "K_A",
                     "alpha_opt", "beta_opt", "t_accl_opt", "alpha_crit",
                     "beta_crit", "t_accl_crit", "softplus_s")
  expect_true(all(thermal_knobs %in% names(s)))
  # ATLS defaults
  expect_identical(s$tcrit_0, 38.0)
  expect_identical(s$k_r1_0, 864.0)
})

test_that("TF24t appends two acclimation ODE states after TF24's states", {
  p <- TF24t_Individual()
  nm <- p$ode_names
  expect_true(all(c("acclim_topt", "acclim_tcrit") %in% nm))
  # the two acclimation states are appended last (inherited indices unchanged)
  expect_identical(tail(nm, 2), c("acclim_topt", "acclim_tcrit"))
  expect_identical(length(nm), length(TF24_Individual()$ode_names) + 2L)
})

test_that("TF24t validates acclimation kinetics", {
  p <- scm_base_parameters("TF24t") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  s <- p$strategies[[1]]; s$alpha_crit <- -1; p$strategies[[1]] <- s
  expect_error(run_scm(p, collect = FALSE, refine_schedule = FALSE),
               "acclimation kinetics")
})

test_that("TF24t runs end-to-end; acclimation sits at its env equilibrium", {
  p <- scm_base_parameters("TF24t") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  res <- run_scm(p, collect = TRUE, refine_schedule = FALSE)
  sp <- res$species  # tidy tibble; state variables are columns
  expect_true(all(c("acclim_topt", "acclim_tcrit") %in% names(sp)))
  expect_true(all(is.finite(sp$acclim_tcrit)))
  expect_true(all(sp$acclim_tcrit >= -1e-8))

  # Constant environment (default leaf_temp = 25) => each state sits at its
  # seeded equilibrium A_eq = (alpha/beta) * softplus(T - T_accl).
  s <- p$strategies[[1]]
  softplus <- function(x, k) log1p(exp(k * x)) / k
  a_crit_eq <- s$alpha_crit / s$beta_crit *
    softplus(25 - s$t_accl_crit, s$softplus_s)
  expect_equal(unique(round(sp$acclim_tcrit, 10)), round(a_crit_eq, 10))
})

test_that("with zero gain the acclimation states stay put (no drift)", {
  p <- scm_base_parameters("TF24t") |>
    add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
  s <- p$strategies[[1]]
  s$alpha_opt <- 0; s$beta_opt <- 0; s$alpha_crit <- 0; s$beta_crit <- 0
  p$strategies[[1]] <- s
  res <- run_scm(p, collect = TRUE, refine_schedule = FALSE)
  # seeded at 0 (alpha=0) and no dynamics -> stays 0
  expect_true(all(abs(res$species$acclim_topt) < 1e-8))
  expect_true(all(abs(res$species$acclim_tcrit) < 1e-8))
})
