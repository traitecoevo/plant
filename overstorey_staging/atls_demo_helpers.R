## Shared leaf- and stand-driving helpers for the TF24t ATLS demo (thermal
## damage, repair, acclimation and avoidance; issue #566).
##
## Kept separate from the .qmd so the smoke test (tests/testthat/test-atls-leaf-demo.R)
## exercises exactly the code the demo runs, and the demo cannot silently rot.
## No plotting here.
##
## Assumes `plant` is loaded (library(plant) or pkgload::load_all()). The layer
## lives on phylloptim's Leaf behind use_thermal_damage_; TF24t is the strategy
## that turns it on. The leaf sections drive a bare Leaf so the four strategy axes
## can be dialled independently.
##
## Rate constants are per day throughout. The `leaf_temp` driver is the AIR
## temperature on the energy-balance path, which TF24t always takes.

## A Leaf with TF24's defaults, as TF24_Strategy::prepare_strategy builds it.
atls_make_leaf <- function() {
  s <- TF24_Strategy(); p <- s$pars; ctrl <- Control()
  Leaf(vcmax_25 = p$vcmax_25, jmax_25 = p$jmax_25, stem_c = p$stem_c,
       stem_P50 = p$stem_P50, root_c = p$root_c, root_P50 = p$root_P50,
       TF24_beta2 = p$TF24_beta2, phi_J = p$phi_J,
       curv_fact_elec_trans = p$curv_fact_elec_trans,
       curv_fact_colim = p$curv_fact_colim,
       GSS_tol_abs = ctrl$GSS_tol_abs,
       vulnerability_curve_ncontrol = ctrl$vulnerability_curve_ncontrol,
       ci_abs_tol = ctrl$ci_abs_tol, ci_niter = ctrl$ci_niter,
       TF24_cost_scale = p$TF24_cost_scale)
}

## The strategy axes as leaf inputs. `atls_traits()` is the default
## ("generalist") set; override any entry through `...`.
##   tolerance   -> topt_offset: one thermostability offset that moves T_opt AND
##                  the impairment's isokinetic temperature together; the acclimated
##                  parts saturate at dTopt_max and dTimpair_max, half-saturation
##                  K_A
##   acclimation -> A: the inducible thermostability state, here held at a
##                  standing value to represent an acclimated leaf
##   repair      -> preventive k_i (rate of functional -> impaired at the
##                  isokinetic temperature k_i_Tiso) and restorative k_rec
##                  (impaired -> functional); k_mat (impaired -> damaged) sets
##                  the lasting floor
##   impairment  -> I, D: the impaired and damaged fractions of leaf area, set
##                  directly on a snapshot leaf; capacity scales by 1 - I - D
## Avoidance is not a leaf trait: it is transpirational cooling on the energy-
## balance path, dialled through the leaf configuration.
atls_traits <- function(...) {
  tr <- list(topt_offset = 0, k_i = 276, k_i_Ea = 316000, k_i_Tiso = 50.5,
             k_rec = 0.25, k_rec_Ea = 0, k_mat = 0.02, k_mat_Ea = 0,
             dTopt_max = 6, dTimpair_max = 6, K_A = 1, A = 0, I = 0, D = 0)
  utils::modifyList(tr, list(...))
}

## Named archetypes, each a single-axis move off the generalist.
atls_strategies <- function() {
  list(
    "Generalist"  = atls_traits(),
    "Tolerant"    = atls_traits(topt_offset = 6),
    "Repairer"    = atls_traits(k_rec = 2.0, k_i = 55.2),
    "Acclimated"  = atls_traits(A = 5)
  )
}

## Hand a trait set to a leaf and turn the layer on.
atls_apply_traits <- function(l, tr = atls_traits()) {
  l$use_thermal_damage_ <- TRUE
  l$topt_offset_  <- tr$topt_offset
  l$k_i_          <- tr$k_i
  l$k_i_Ea_       <- tr$k_i_Ea
  l$k_i_Tiso_     <- tr$k_i_Tiso
  l$k_rec_        <- tr$k_rec
  l$k_rec_Ea_     <- tr$k_rec_Ea
  l$k_mat_        <- tr$k_mat
  l$k_mat_Ea_     <- tr$k_mat_Ea
  l$dTopt_max_    <- tr$dTopt_max
  l$dTimpair_max_ <- tr$dTimpair_max
  l$K_A_          <- tr$K_A
  l$A_            <- tr$A
  l$I_            <- tr$I
  l$D_            <- tr$D
  invisible(l)
}

## A well-watered, well-rooted, moderate-conductance operating point, so the
## profit optimum opens the stomata (a shut leaf makes the contrasts
## degenerate). A better-coupled "avoider" uses a smaller leaf dimension d,
## higher conductance or more wind.
atls_leaf_config <- function(...) {
  cfg <- list(root_carbon_per_leaf_area = 20.0, psi_soil = 0.3,
              leaf_specific_conductance_max = 5e-4,
              ca = 40, atm_o2_kpa = 21, atm_kpa = 101.3, d = 0.05,
              wind_speed = 2.0)
  utils::modifyList(cfg, list(...))
}

## Seat a leaf at one environment cell with a trait set. `Tenv` is the
## prescribed leaf temperature off the energy balance (pm = FALSE) and the air
## temperature on it (pm = TRUE).
atls_set_physiology <- function(l, PAR, Tenv, VPD, tr = atls_traits(),
                                pm = FALSE, cfg = atls_leaf_config()) {
  atls_apply_traits(l, tr)
  l$use_energy_balance_ <- isTRUE(pm)
  l$d_ <- cfg$d
  l$wind_speed_ <- cfg$wind_speed
  # The intercepted near-infrared in daylight: as much energy as the PAR band.
  l$NIR_int_ <- PAR / 4.57
  # The leaf takes resistances, so the root architecture model runs here, with
  # the constants TF24_Strategy uses.
  l$set_physiology(
    root_network = phylloptim::root_network_from_carbon(
      cfg$root_carbon_per_leaf_area, soil_depth = 1,
      beta_R_H = 3.4e2, beta_R_V = 9.4e3),
    PPFD = PAR, psi_soil = cfg$psi_soil, soil_depth = 1,
    leaf_specific_conductance_max = cfg$leaf_specific_conductance_max,
    atm_vpd = VPD, ca = cfg$ca,
    leaf_temp = Tenv, atm_o2_kpa = cfg$atm_o2_kpa, atm_kpa = cfg$atm_kpa)
  invisible(l)
}

## The impairment rate k_imp(Tleaf), per day, for a trait set: k_i at the
## isokinetic temperature k_i_Tiso (50.5 deg C at the defaults) plus the
## tolerance and acclimation shifts, and an Arrhenius rise of tenfold about every
## 6 K. The leaf temperature is prescribed, so this is the mechanism with no
## optimiser in the loop.
atls_k_imp_curve <- function(Tleaf_seq, tr = atls_traits(), PAR = 1500, VPD = 2,
                             cfg = atls_leaf_config()) {
  vapply(Tleaf_seq, function(T) {
    l <- atls_make_leaf()
    atls_set_physiology(l, PAR, T, VPD, tr, pm = FALSE, cfg)
    l$k_imp_eff_
  }, numeric(1))
}

## The standing impaired fraction the I ODE relaxes toward at a held leaf
## temperature with D = 0,
##   I* = k_imp / (k_imp + k_rec_eff + k_mat),
## from the leaf's own rates (per day; the day-to-year factor cancels).
atls_I_star_curve <- function(Tleaf_seq, tr = atls_traits(), PAR = 1500,
                              VPD = 2, cfg = atls_leaf_config()) {
  vapply(Tleaf_seq, function(T) {
    l <- atls_make_leaf()
    atls_set_physiology(l, PAR, T, VPD, tr, pm = FALSE, cfg)
    l$k_imp_eff_ / (l$k_imp_eff_ + l$k_rec_eff_ + tr$k_mat)
  }, numeric(1))
}

## Solve the operating point at one environment cell; a 1-row data.frame of
## leaf outputs. Tleaf is the leaf's own operating temperature: the prescribed
## one off the energy balance, the solved one on it.
atls_solve_cell <- function(PAR, Tenv, VPD, tr = atls_traits(), pm = FALSE,
                            cfg = atls_leaf_config(), label = NA_character_) {
  l <- atls_make_leaf()
  atls_set_physiology(l, PAR, Tenv, VPD, tr, pm, cfg)
  l$find_root_collar_psi()
  k_imp <- l$k_imp_eff_
  data.frame(strategy = label, PAR = PAR, Tenv = Tenv, VPD = VPD, pm = pm,
             Tleaf = l$Tleaf_, k_imp = k_imp,
             I_star = k_imp / (k_imp + l$k_rec_eff_ + tr$k_mat),
             A = l$assim_colimited_, gs = l$stom_cond_CO2_,
             E = l$transpiration_, profit = l$profit_)
}

## Solve a Tenv gradient for every named strategy (off the energy balance by
## default); a long data.frame with a `strategy` column.
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
## A single-species patch under a CONSTANT climate (the `leaf_temp` driver as the
## midday air temperature; other drivers at their defaults), run for a whole
## patch lifetime. Each run returns the net reproduction ratio R0 and lifetime
## offspring production. Runs take tens of seconds, so callers keep grids small
## and the .qmd caches them.

## One SCM run. `type` is "TF24t", or "TF24" for the energy-balance-only
## comparator (pars$use_energy_balance set to 1 to match TF24t's forced
## balance). `mutate` edits the resident strategy: the thermal axes are TF24t
## strategy fields (topt_offset for tolerance; k_i and k_rec for repair; alpha
## for the acclimation kinetics), so a strategy here is a heritable trait set.
atls_scm_fitness <- function(air_temp, type = "TF24t", mutate = identity,
                             lma = 0.0825, birth_rate = 20, PPFD = NULL) {
  p <- scm_base_parameters(type) |>
    add_strategies(trait_matrix(lma, "lma"), birth_rate = birth_rate)
  s <- p$strategies[[1]]
  if (identical(type, "TF24")) s$pars$use_energy_balance <- 1
  s <- mutate(s)
  p$strategies[[1]] <- s
  env <- Environment(type)
  env$extrinsic_drivers_set_constant("leaf_temp", air_temp)
  if (!is.null(PPFD)) env$extrinsic_drivers_set_constant("SWdown", PPFD / (0.5 * 4.57))
  res <- run_scm(p, env = env, collect = FALSE, refine_schedule = FALSE)
  data.frame(type = type, air_temp = air_temp,
             R0 = res$net_reproduction_ratios,
             offspring = res$offspring_production)
}

## Heritable archetypes for the SCM comparisons, each a single-axis move off the
## generalist. Acclimation is a change of kinetics, since A is a state.
atls_scm_strategies <- function() {
  list(
    "Generalist" = identity,
    "Tolerant"   = function(s) { s$topt_offset <- 6; s },
    "Repairer"   = function(s) { s$k_rec <- 2.0; s$k_i <- 55.2; s },
    "Acclimator" = function(s) { s$alpha <- 0.06; s })
}

## One worker per core, one held back for the session.
atls_scm_default_workers <- function() {
  n <- parallel::detectCores()
  if (is.na(n)) n <- 1L
  max(1L, n - 1L)
}

## Run a list of zero-argument jobs, each returning a 1-row data.frame, and rbind
## them in order. Forked (parallel::mclapply), so the workers inherit whatever
## plant the session loaded -- a pkgload::load_all build included, which a
## PSOCK cluster would not see. Sequential on Windows or with one worker. SCM
## runs are deterministic and mclapply keeps order, so both give the same frame.
atls_scm_run_jobs <- function(jobs, workers = atls_scm_default_workers()) {
  workers <- max(1L, min(workers, length(jobs)))
  use_fork <- workers > 1L && .Platform$OS.type != "windows"
  rows <- if (use_fork) {
    parallel::mclapply(jobs, function(f) f(), mc.cores = workers,
                       mc.preschedule = FALSE)
  } else {
    lapply(jobs, function(f) f())
  }
  failed <- vapply(rows, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop("SCM job failed: ", as.character(rows[[which(failed)[1]]]))
  }
  do.call(rbind, rows)
}

## --- Fitness: viability, not R0 at a fixed birth rate -------------------------
## R0 at a shared, non-trivial birth rate is NOT a valid cross-strategy measure:
## it conflates per-capita productivity with density-dependent recruitment. A
## less fecund resident shades the patch less, so at a shared seed rain each of
## its seeds does better, inflating its R0 although it is the weaker plant.
##
## The density-free measure is viability: does a rare resident more than replace
## itself? R0 is evaluated in the low-density limit (a small birth rate). R0 >= 1
## means the population grows from rare; the warmest climate with R0 >= 1 is the
## persistence boundary, which is comparable across strategies because density
## effects vanish as density goes to zero. Read the boundary, never the R0
## magnitude away from it.

## The small birth rate standing in for the low-density limit.
atls_low_density_birth_rate <- 1

atls_or <- function(a, b) if (is.null(a)) b else a

## Low-density R0 across a climate gradient for named scenarios, each a
## list(type =, mutate =); the `scenario` column labels them.
atls_scm_viability <- function(air_temp_seq, scenarios,
                               birth_rate = atls_low_density_birth_rate,
                               workers = atls_scm_default_workers()) {
  grid <- expand.grid(lt = air_temp_seq, nm = names(scenarios),
                      stringsAsFactors = FALSE)
  jobs <- lapply(seq_len(nrow(grid)), function(i) {
    nm <- grid$nm[i]; lt <- grid$lt[i]; sc <- scenarios[[nm]]
    function() {
      r <- atls_scm_fitness(lt, type = sc$type,
                            mutate = atls_or(sc$mutate, identity),
                            birth_rate = birth_rate)
      r$scenario <- nm
      r
    }
  })
  atls_scm_run_jobs(jobs, workers)
}

## The persistence boundary per scenario: the warmest air temperature with
## low-density R0 >= 1, by log-linear interpolation across the sweep. The
## collapse is steep, so the boundary is good to about half a grid step. NA where
## the scenario never crosses 1 inside the swept range.
atls_persistence_boundary <- function(viab) {
  do.call(rbind, lapply(split(viab, viab$scenario), function(d) {
    d <- d[order(d$air_temp), ]
    lr <- log(pmax(d$R0, 1e-300))
    above <- lr >= 0
    b <- NA_real_
    if (any(above) && any(!above)) {
      i <- which(above[-length(above)] & !above[-1])[1]
      if (!is.na(i)) {
        b <- d$air_temp[i] +
          (0 - lr[i]) / (lr[i + 1] - lr[i]) * (d$air_temp[i + 1] - d$air_temp[i])
      }
    }
    data.frame(scenario = d$scenario[1], boundary = b)
  }))
}
