## Shared leaf-driving helpers for the TF24t ATLS leaf demo (thermal
## damage / repair / acclimation / avoidance; issue #566).
##
## Kept separate from the .qmd so the smoke test (tests/testthat/test-atls-leaf-demo.R)
## exercises exactly the same leaf-driving code as the vignette, guaranteeing the
## demo cannot silently rot. No plotting here -- pure model driving.
##
## Assumes the `plant` package is loaded (library(plant) or devtools::load_all()).
## The ATLS layer lives on the base Leaf (use_thermal_damage_ + the thermal
## traits); TF24t is the strategy that turns it on. This demo drives the Leaf
## directly so the four strategy axes can be dialled independently.

## PM leaf-energy-balance physical constants (mirror inst/include/plant/leaf_model.h),
## used only to recover the operating-point leaf temperature Tleaf in R.
atls_lambda <- 2.45e6   # latent heat of vaporisation, J kg^-1
atls_rho_cp <- 1200.0   # volumetric heat capacity of air, J m^-3 K^-1

## Build a Leaf with TF24 defaults, matching TF24_Strategy::prepare_strategy.
atls_make_leaf <- function() {
  s <- TF24_Strategy(); p <- s$pars; ctrl <- Control()
  root_c <- 2.680147; root_b <- 3.898245
  root_psi_crit <- root_b * (log(1 / 0.05))^(1 / root_c)
  Leaf(vcmax_25 = p$vcmax_25, jmax_25 = p$jmax_25, c = p$c, b = p$b,
       psi_crit = p$psi_crit, root_c = root_c, root_b = root_b,
       root_psi_crit = root_psi_crit, beta2 = p$beta2, a = p$a,
       curv_fact_elec_trans = p$curv_fact_elec_trans,
       curv_fact_colim = p$curv_fact_colim,
       GSS_tol_abs = ctrl$GSS_tol_abs,
       vulnerability_curve_ncontrol = ctrl$vulnerability_curve_ncontrol,
       ci_abs_tol = ctrl$ci_abs_tol, ci_niter = ctrl$ci_niter,
       g1_TF24 = p$g1_TF24, beta_R_H = 3.4e2, beta_R_V = 9.4e3)
}

## The ATLS merged-revision strategy axes as leaf traits. `atls_traits()` returns
## the default ("generalist") trait set; override any knob via ....
##   tolerance   -> topt_offset (one thermostability offset shifting T_opt AND the
##                  damage onset together), with ceiling dTopt_max / half-sat K_A
##   acclimation -> A (inducible thermostability state; here set to a standing
##                  value to represent a fully-acclimated leaf)
##   repair      -> preventive k_i (leak phi_d->I_r) and restorative k_rec
##                  (resynthesis of I_r); k_mat sets the permanent floor
##   damage      -> I_r / I_p (the lasting-damage pools; a leaf snapshot sets them
##                  directly to show the (1 - I_r - I_p) capacity discount)
##   (avoidance is not a leaf trait -- it is transpirational cooling on the PM
##    path; dial it via the leaf config: conductance / dimension d / wind.)
atls_traits <- function(...) {
  tr <- list(topt_offset = 0, k_i = 0.5, k_rec = 0.25, k_mat = 0.02,
             m_rep = 0.4, t_rep_cut = 45, dTopt_max = 6, K_A = 1,
             A = 0, I_r = 0, I_p = 0)
  utils::modifyList(tr, list(...))
}

## Named archetypes spanning the strategy space (each a single-axis move off the
## generalist so the plots read cleanly).
atls_strategies <- function() {
  list(
    "Generalist"  = atls_traits(),
    "Tolerant"    = atls_traits(topt_offset = 6),
    "Repairer"    = atls_traits(k_rec = 2.0, k_i = 0.1),
    "Acclimated"  = atls_traits(A = 5)
  )
}

## Apply a trait set to a leaf and turn the damage layer on.
atls_apply_traits <- function(l, tr = atls_traits()) {
  l$use_thermal_damage_ <- TRUE
  l$topt_offset_ <- tr$topt_offset
  l$k_i_         <- tr$k_i
  l$k_rec_       <- tr$k_rec
  l$k_mat_       <- tr$k_mat
  l$m_rep_       <- tr$m_rep
  l$t_rep_cut_   <- tr$t_rep_cut
  l$dTopt_max_   <- tr$dTopt_max
  l$K_A_         <- tr$K_A
  l$A_           <- tr$A
  l$I_r_         <- tr$I_r
  l$I_p_         <- tr$I_p
  invisible(l)
}

## Well-watered, well-rooted, moderate-conductance operating point (so the profit
## optimum opens the stomata; a shut-down leaf makes contrasts degenerate).
## Override via ... -- e.g. a better-coupled "avoider" uses larger conductance /
## smaller leaf dimension d / higher wind.
atls_leaf_config <- function(...) {
  cfg <- list(area_leaf = 1.0, mass_root_prop = 20.0, psi_soil = 0.3,
              leaf_specific_conductance_max = 5e-3,
              sapwood_volume_per_leaf_area = 1e-3,
              ca = 40, atm_o2_kpa = 21, atm_kpa = 101.3, d = 0.05, wind_speed = 2.0)
  utils::modifyList(cfg, list(...))
}

## Configure physiology for one environment cell, with a trait set. `Tenv` is the
## prescribed leaf temperature on the Fick path (pm = FALSE) or the air
## temperature on the PM path (pm = TRUE).
atls_set_physiology <- function(l, PAR, Tenv, VPD, tr = atls_traits(),
                                pm = FALSE, cfg = atls_leaf_config()) {
  atls_apply_traits(l, tr)
  l$use_energy_balance_ <- isTRUE(pm)
  l$d_ <- cfg$d
  l$wind_speed_ <- cfg$wind_speed
  l$set_physiology(
    area_leaf = cfg$area_leaf, mass_root_prop = cfg$mass_root_prop, rho = 608,
    a_bio = 0.0245, PPFD = PAR, psi_soil = cfg$psi_soil, soil_depth = 1,
    leaf_specific_conductance_max = cfg$leaf_specific_conductance_max,
    atm_vpd = VPD, ca = cfg$ca,
    sapwood_volume_per_leaf_area = cfg$sapwood_volume_per_leaf_area,
    leaf_temp = Tenv, atm_o2_kpa = cfg$atm_o2_kpa, atm_kpa = cfg$atm_kpa)
  invisible(l)
}

## Operating-point leaf temperature (deg C). On the PM path, recovered from the
## exposed energy-balance fields; on the Fick path Tleaf == Tenv by definition.
atls_leaf_temp <- function(l, Tenv) {
  if (isTRUE(l$use_energy_balance_)) {
    l$Tair_ + (l$Rn_ - atls_lambda * l$transpiration_) * l$ra_ / atls_rho_cp
  } else {
    Tenv
  }
}

## The reversible deactivation curve phi_d(Tleaf) for a trait set: the
## transiently-unfolded fraction, K/(1+K), from the SAME Medlyn K(T) inside the
## peaked-Arrhenius curve. It is the substrate feeding the irreversible I_r leak;
## its half-unfolding (onset) is emergent at T_K = H_d/d_S ~ 34.5 C. Prescribed
## Tleaf (Fick path) so it is the pure mechanism, computed in set_physiology and
## read straight back.
atls_phi_d_curve <- function(Tleaf_seq, tr = atls_traits(), PAR = 1500, VPD = 2,
                             cfg = atls_leaf_config()) {
  vapply(Tleaf_seq, function(T) {
    l <- atls_make_leaf()
    atls_set_physiology(l, PAR, T, VPD, tr, pm = FALSE, cfg)
    l$phi_d_
  }, numeric(1))
}

## Illustrative quasi-steady recoverable-damage fraction I_r*(Tleaf): the balance
## the I_r ODE relaxes toward at a held leaf temperature, with I_p = 0,
##   I_r* = k_i * phi_d / (k_i * phi_d + k_rec_eff + k_mat),
## using the leaf's own phi_d and gated k_rec_eff (day^-1; the DAYS_PER_YEAR
## factor cancels). Shows how tolerance (lower phi_d) and repair (higher k_rec)
## reduce standing recoverable damage.
atls_Ir_star_curve <- function(Tleaf_seq, tr = atls_traits(), PAR = 1500,
                               VPD = 2, cfg = atls_leaf_config()) {
  vapply(Tleaf_seq, function(T) {
    l <- atls_make_leaf()
    atls_set_physiology(l, PAR, T, VPD, tr, pm = FALSE, cfg)
    influx <- tr$k_i * l$phi_d_
    influx / (influx + l$k_rec_eff_ + tr$k_mat)
  }, numeric(1))
}

## Solve the optimal operating point at one environment cell for a trait set;
## returns a 1-row data.frame of leaf outputs (N, assimilation, Tleaf, ...).
## After find_root_collar_psi, N_ / assim_colimited_ reflect the operating point
## (Fick: at the prescribed Tleaf; PM: at the cooled operating-point Tleaf).
atls_solve_cell <- function(PAR, Tenv, VPD, tr = atls_traits(), pm = FALSE,
                            cfg = atls_leaf_config(), label = NA_character_) {
  l <- atls_make_leaf()
  atls_set_physiology(l, PAR, Tenv, VPD, tr, pm, cfg)
  l$find_root_collar_psi()
  influx <- tr$k_i * l$phi_d_
  data.frame(strategy = label, PAR = PAR, Tenv = Tenv, VPD = VPD, pm = pm,
             Tleaf = atls_leaf_temp(l, Tenv), phi_d = l$phi_d_,
             Ir_star = influx / (influx + l$k_rec_eff_ + tr$k_mat),
             A = l$assim_colimited_, gs = l$stom_cond_CO2_,
             E = l$transpiration_, profit = l$profit_)
}

## Solve a Tenv gradient for every named strategy (Fick path by default), long
## data.frame with a `strategy` column. Used for the strategy-comparison plots.
atls_solve_strategies <- function(Tenv_seq, strategies = atls_strategies(),
                                  PAR = 1500, VPD = 2, pm = FALSE,
                                  cfg = atls_leaf_config()) {
  rows <- list(); k <- 0L
  for (nm in names(strategies)) {
    for (Tenv in Tenv_seq) {
      k <- k + 1L
      rows[[k]] <- atls_solve_cell(PAR, Tenv, VPD, strategies[[nm]], pm, cfg,
                                   label = nm)
    }
  }
  do.call(rbind, rows)
}

## --- Community scale (SCM) --------------------------------------------------
## The leaf sections above are a single isolated leaf. These drive the full
## size-structured community model: a single-species patch under a *constant*
## climate (the air_temp driver = the midday air temperature; other drivers at
## their defaults) run to demographic equilibrium, returning community-fitness
## scalars -- the net reproduction ratio R0 (a resident persists iff R0 >= 1) and
## lifetime offspring production. SCM runs are ~30 s each, so callers keep the
## grids small; the .qmd caches its chunks.

## One SCM run. `type` is "TF24t" (PM + ATLS) or "TF24" (a PM-only comparator,
## with use_energy_balance forced on to match TF24t's forced PM). `mutate` tweaks
## the resident strategy -- the thermal axes are TF24t *strategy* fields
## (topt_offset for constitutive thermostability; k_i / k_rec for repair; alpha
## for the acclimation kinetics), so an SCM strategy is a heritable trait set,
## unlike the fixed-A leaf snapshots in atls_strategies().
atls_scm_fitness <- function(air_temp, type = "TF24t", mutate = identity,
                             lma = 0.0825, birth_rate = 20, PPFD = NULL) {
  p <- scm_base_parameters(type) |>
    add_strategies(trait_matrix(lma, "lma"), birth_rate = birth_rate)
  s <- p$strategies[[1]]
  if (identical(type, "TF24")) s$pars$use_energy_balance <- 1  # PM-only comparator
  s <- mutate(s)
  p$strategies[[1]] <- s
  env <- Environment(type)
  env$extrinsic_drivers_set_constant("air_temp", air_temp)
  if (!is.null(PPFD)) env$extrinsic_drivers_set_constant("PPFD", PPFD)
  res <- run_scm(p, env = env, collect = FALSE, refine_schedule = FALSE)
  data.frame(type = type, air_temp = air_temp,
             R0 = res$net_reproduction_ratios,
             offspring = res$offspring_production)
}

## Heritable thermal-strategy archetypes for the SCM tournament (each a
## single-axis move off the generalist; acclimation is a *kinetics* change since
## it is a dynamic state, not a fixed trait).
atls_scm_strategies <- function() {
  list(
    "Generalist" = identity,
    "Tolerant"   = function(s) { s$topt_offset <- 6; s },
    "Repairer"   = function(s) { s$k_rec <- 2.0; s$k_i <- 0.1; s },
    "Acclimator" = function(s) { s$alpha <- 0.06; s })
}

## R0 across a climate gradient for TF24t vs the PM-only TF24 comparator.
atls_scm_climate <- function(air_temp_seq, types = c("TF24", "TF24t")) {
  rows <- list(); k <- 0L
  for (ty in types) for (lt in air_temp_seq) {
    k <- k + 1L
    rows[[k]] <- atls_scm_fitness(lt, type = ty)
  }
  do.call(rbind, rows)
}

## R0 across a climate gradient for every heritable strategy archetype (TF24t).
atls_scm_tournament <- function(air_temp_seq, strategies = atls_scm_strategies()) {
  rows <- list(); k <- 0L
  for (nm in names(strategies)) for (lt in air_temp_seq) {
    k <- k + 1L
    r <- atls_scm_fitness(lt, type = "TF24t", mutate = strategies[[nm]])
    r$strategy <- nm
    rows[[k]] <- r
  }
  do.call(rbind, rows)
}
