# TF24f stand-metric gradients (#472 scope B, build-order step 1). TF24f is the
# fast-acclimation TF24 variant whose optimal root-collar potential is a 6th ODE state
# tracked by gradient ascent (no per-step golden-section optimiser). TF24f -- not TF24
# -- is the right target for the
# stand CENSUS gradients: its leaf evaluation at the tracked collar is analytic / IFT-able,
# so the census number density's growth-rate-gradient term needs no curvature harvest. The
# R0 GATE below is the first deliverable: a double-precision census reconstruction that
# proves the collar-state replay is faithful (the prerequisite §7 flags) before the
# reverse-mode tape (R1) is built.

# Harvest a run-with-cache TF24f SCM into the pieces the census replay consumes (the
# TF24f mirror of tf24_harvest). Unlike TF24, no per-stage leaf-opt harvest is needed --
# the replay evaluates the analytic leaf at the tracked collar live -- so this returns
# only the frozen schedule + env, the cohort birth steps, the parameter vector and the
# acclimation knobs (k_acclim / use_ad_gradient). Requires shading_model = "crown-centre".
tf24f_harvest <- function(scm, species = 1L, birth_rate = NULL) {
  types <- extract_RcppR6_template_types(scm$parameters, "Parameters")
  if (!identical(types[[1]], "TF24f")) {
    stop("TF24f census gradients are implemented for the TF24f strategy only")
  }
  if (species < 1L || species > length(scm$patch$species)) {
    stop("species index out of range: stand has ", length(scm$patch$species),
         " species")
  }
  patch <- scm$patch
  sh <- patch$step_history
  eh <- patch$environment_history
  if (length(eh) < 1L) {
    stop("No resident schedule cached: run the SCM with control(save_RK45_cache = TRUE)")
  }
  sp    <- patch$species[[species]]
  nt    <- sp$node_times
  pdens <- sp$patch_densities
  ppsab <- sp$pr_patch_survival_at_birth
  strat <- scm$parameters$strategies[[species]]
  pp    <- unlist(strat$pars)

  if (is.null(birth_rate)) {
    birth_rate <- scm$offspring_production[[species]] / scm$net_reproduction_ratios[[species]]
  }
  # Cohort birth steps (introductions land on step times).
  birth_step <- vapply(nt, function(t) which.min(abs(sh - t)) - 1L, integer(1))
  N <- length(eh)

  # offspring_production weighting (the TF24f mirror of ff16_harvest / tf24_harvest): the
  # node-spacing trapezoid weights so offspring_production == sum_i tw_i * offspring_i,
  # and per-RK-stage patch survival at the exact Cash-Karp stage times. Only the offspring
  # tape uses these (census / resident do not).
  tcoef <- numeric(length(nt)); x <- nt; n <- length(x)
  tcoef[1] <- 0.5 * (x[2] - x[1]); tcoef[n] <- 0.5 * (x[n] - x[n - 1])
  if (n > 2) tcoef[2:(n - 1)] <- 0.5 * (x[3:n] - x[1:(n - 2)])
  tw <- tcoef * pdens * pp[["S_D"]] * birth_rate
  ah <- c(0, 0.2, 0.3, 0.6, 1.0, 0.875); hN <- diff(sh)
  ppsurv <- matrix(0, N, 6)
  for (k in seq_len(N)) for (s in 1:6) ppsurv[k, s] <- patch$pr_survival(sh[k] + ah[s] * hN[k])

  # Boundary new_node (height + competition effect) per RK stage: the trapezium tail
  # term Species::compute_competition adds beyond `nodes`, needed by the COUPLED
  # resident replay's per-stage canopy reconstruction at ground level (the TF24f mirror
  # of ff16_harvest's nn_h / nn_c). Empty on older caches.
  nn_h <- patch$stand_newnode_height_stage_history
  nn_c <- patch$stand_newnode_competition_stage_history
  patch_area <- scm$parameters$patch_area

  list(pp = pp, eh = eh, sh = sh, birth_step = birth_step, birth_rate = birth_rate,
       k_acclim = strat$k_acclim, use_ad_gradient = strat$use_ad_gradient, nt = nt,
       tw = tw, ppsurv = ppsurv, ppsab = ppsab,
       nn_h = nn_h, nn_c = nn_c, patch_area = patch_area)
}

# R0 gate (internal): double-precision census reconstruction of a TF24f resident stand.
# Re-evolves every cohort's {5 demog, tracked collar, log_density} over the frozen
# schedule and returns the reconstructed per-cohort heights / collar / log-densities and
# the census metric values. Used to confirm the collar-state replay reproduces the SCM's
# stored stand (heights / log_densities / opt_root_psi_state) and that the LAI reduction
# matches compute_competition(0), before the reverse-mode tape (R1) is added. Not yet a
# public gradient entry point -- it returns the recon, not d(metric)/d(theta).
# Strategy guard for the native census entries (which take the live SCM via an
# RcppR6<SCM<TF24f,...>> cast). Mirrors the check tf24f_harvest used to provide, so a
# non-TF24f SCM still gets a clear message rather than a raw external-pointer cast error.
tf24f_require_strategy <- function(scm) {
  types <- extract_RcppR6_template_types(scm$parameters, "Parameters")
  if (!identical(types[[1]], "TF24f"))
    stop("TF24f census gradients are implemented for the TF24f strategy only")
}

tf24f_census_recon <- function(scm, metrics = c("LAI", "biomass", "size_moment"),
                               species = 1L, birth_rate = NULL) {
  tf24f_require_strategy(scm)
  # The reconstruction's growth-rate gradient g' must match whatever the resident
  # SCM used (Node::growth_rate_gradient): the exact-AD path when the run set
  # control(node_gradient_exact_ad = TRUE), else the backward finite difference.
  exact_ad <- isTRUE(scm$parameters$strategies[[species]]$control$node_gradient_exact_ad)
  # Fully native: env + birth steps from the live Patch (no Rcpp::as<> round-trip; the
  # whole harvest is native, so the R-side tf24f_harvest -- whose ppsurv/tw loop is dead
  # work for census -- is skipped). birth_rate < 0 recovers the rate natively.
  strat <- scm$parameters$strategies[[species]]
  tf24f_census_recon_native(scm, unlist(strat$pars), as.integer(species - 1L),
                            if (is.null(birth_rate)) -1 else birth_rate,
                            strat$k_acclim, strat$use_ad_gradient, metrics, exact_ad)
}

# TF24f frozen census trait gradient (#472 scope B, build-order step 2 -- R1 GATE).
# d(census metric)/d(theta) for the FROZEN (rare-mutant / invasion) resident light, by a
# central finite difference over the R0 census reconstruction: for each trait, perturb the
# parameter vector by +/- a relative step, re-run the double-precision collar-state replay
# against the SAME frozen resident environment, and difference the metric. The resident
# schedule + per-RK-stage light are held fixed (the invasion gradient), so this is the
# faithful finite-difference reference the reverse-mode AD tape (the actual R1) must
# reproduce -- and a usable prototype gradient in its own right. Returns a metrics x traits
# Jacobian and the reconstructed metric values. (The AD tape will replace the per-trait
# replays with one reverse sweep; this gate fixes the target it must hit.)
tf24f_census_gradient_fd <- function(scm, metrics = c("LAI", "biomass", "size_moment"),
                                     traits = NULL, species = 1L, birth_rate = NULL,
                                     rel_step = 1e-5) {
  if (is.null(traits)) traits <- tf24_default_traits()
  h <- tf24f_harvest(scm, species, birth_rate)
  exact_ad <- isTRUE(scm$parameters$strategies[[species]]$control$node_gradient_exact_ad)
  # Native recon so the FD reference differentiates the SAME function the native AD
  # tape does (faithful crown-sampled light, no Rcpp::as<> round-trip). birth steps are
  # computed natively inside; h$birth_rate is the concrete recovered rate.
  recon <- function(pp) tf24f_census_recon_native(scm, pp, as.integer(species - 1L),
    h$birth_rate, h$k_acclim, h$use_ad_gradient, metrics, exact_ad)$values

  values <- recon(h$pp)
  jac <- matrix(0, length(metrics), length(traits),
                dimnames = list(metrics, traits))
  for (j in seq_along(traits)) {
    tr <- traits[j]
    if (!tr %in% names(h$pp)) stop("unknown TF24f trait: ", tr)
    d <- rel_step * max(abs(h$pp[[tr]]), 1e-8)
    pp_p <- h$pp; pp_p[[tr]] <- pp_p[[tr]] + d
    pp_m <- h$pp; pp_m[[tr]] <- pp_m[[tr]] - d
    jac[, j] <- (recon(pp_p) - recon(pp_m)) / (2 * d)
  }
  list(jacobian = jac, values = values)
}

# TF24f frozen census trait gradient by reverse-mode AD (#472 scope B, build-order
# step 2 -- the R1 tape, the refine of tf24f_census_gradient_fd). One reverse sweep over
# the 7-state replay {5 demog, tracked collar, log_density} per metric, replacing the
# per-trait finite difference over the census reconstruction. The hard ingredient (the
# tracked collar is a strongly theta-dependent STATE that lags the optimum, so the
# envelope theorem does not zero it) is handled by a CURVATURE harvest: the collar is
# carried as a taped state with rate k_acclim * dprofit_dpsi, linearised with the second
# derivatives d2profit/dpsi2, d2profit/dpsi dh, d2profit/dpsi dtheta_k (all FD-harvested
# in the double discovery pass). The census g' = d(height_dt)/d(height) reproduces the
# SCM's backward-FD scheme by harvesting a second operating point at h - GEPS. Resident
# light is FROZEN (the rare-mutant / invasion gradient). Validated to ~1% (the recon
# noise floor) against tf24f_census_gradient_fd; see test-tf24f-census-gradient.R.
tf24f_census_gradient_ad <- function(scm, metrics = c("LAI", "biomass", "size_moment"),
                                     traits = NULL, species = 1L, birth_rate = NULL,
                                     trait_rel_step = 1e-5) {
  tf24f_require_strategy(scm)
  if (is.null(traits)) traits <- tf24_default_traits()
  # Fully native: env + birth steps from the live Patch; tf24f_harvest skipped.
  strat <- scm$parameters$strategies[[species]]
  tf24f_census_gradient_ad_native(scm, unlist(strat$pars), as.integer(species - 1L),
                                  if (is.null(birth_rate)) -1 else birth_rate,
                                  strat$k_acclim, strat$use_ad_gradient, traits, metrics,
                                  trait_rel_step)
}

# TF24f individual grow-to-size trait gradient by reverse-mode AD (#472 scope B,
# build-order step 4 -- the AD refine of tf24f_grow_individual_to_size_gradient_fd). A
# single TF24f plant grown in a FIXED environment to target size(s); returns d(t*)/d(theta)
# and the TOTAL d(state at t*)/d(theta). No resident feedback / canopy / density, so the
# only machinery beyond the frozen-schedule replay is the stopping-time IFT. The tracked
# collar (opt_root_psi_state) is carried as the 6th replayed state with the SAME curvature-
# linearised gradient-ascent rate the census tape uses; it starts at the individual's birth
# value (theta-independent), so only the seedling height h0 carries an initial-condition
# derivative. Pass 1 (R) harvests the live grow's adaptive Cash-Karp schedule via
# grow_individual_bracket; pass 2 (C++) replays + sweeps. Validated against
# tf24f_grow_individual_to_size_gradient_fd; see test-tf24f-individual-gradient.R.
tf24f_grow_individual_to_size_gradient_ad <- function(individual, sizes, size_name, env,
                                                      traits = NULL, time_max = Inf,
                                                      warn = TRUE, trait_rel_step = 1e-5) {
  if (!grepl("^TF24f", individual$strategy_name))
    stop("tf24f_grow_individual_to_size_gradient_ad is for the TF24f strategy only")
  if (is.unsorted(sizes) || length(sizes) == 0L)
    stop("sizes must be non-empty and sorted")
  sidx <- match(size_name, individual$ode_names)
  if (is.na(sidx))
    stop("size_name must be one of the ODE state names: ",
         paste(individual$ode_names, collapse = ", "))
  if (is.null(traits)) traits <- tf24_default_traits()

  # Pass 1: harvest the adaptive step schedule (the frozen schedule the replay reproduces).
  brk <- grow_individual_bracket(individual, sizes, size_name, env, time_max, warn)
  y0  <- individual$ode_state
  s   <- individual$strategy
  pp  <- unlist(s$pars)
  shading <- s$control$shading_model
  if (is.null(shading) || !nzchar(shading)) shading <- "mean-light"  # TF24's default

  res <- tf24f_grow_to_size_gradient_impl(pp, env, y0, brk$time, as.numeric(sizes),
                                          as.integer(sidx - 1L), traits, s$k_acclim,
                                          s$use_ad_gradient, shading, s$control$GSS_tol_abs,
                                          trait_rel_step)
  res$sizes <- sizes
  res
}

# TF24f individual grow-to-size trait gradient (#472 scope B, the "individuals" surface
# -- prototype). d(t*)/d(theta) and d(state at t*)/d(theta) for a single plant grown in a
# FIXED environment to target size(s), by central finite difference over
# grow_individual_to_size. No resident feedback (the env is given) and no canopy/density,
# so this is the lightest gradient surface; for TF24f the tracked collar is re-evolved
# inside each grow (it is one of the ODE states), so the FD captures the collar's response
# automatically -- which is exactly why an exact AD version is the heavier follow-up (the
# tracked collar is strongly theta-dependent). The
# FD here is the prototype + the reference that AD version must reproduce. Traits are
# perturbed on the (post-hyperpar) strategy parameters directly, matching the census FD.
tf24f_grow_individual_to_size_gradient_fd <- function(individual, sizes, size_name, env,
                                                      traits = NULL, time_max = Inf,
                                                      warn = FALSE, rel_step = 1e-5) {
  if (!grepl("^TF24f", individual$strategy_name))
    stop("tf24f_grow_individual_to_size_gradient_fd is for the TF24f strategy only")
  if (is.null(traits)) traits <- tf24_default_traits()

  grow <- function(ind) {
    r <- grow_individual_to_size(ind, sizes, size_name, env, time_max, warn)
    list(time = r$time, state = r$state)
  }
  perturb <- function(tr, delta) {
    s <- individual$strategy
    pars <- s$pars
    if (!tr %in% names(pars)) stop("unknown TF24f trait: ", tr)
    pars[[tr]] <- pars[[tr]] + delta
    s$pars <- pars
    grow(TF24f_Individual(s))
  }

  base <- grow(individual)
  nS <- length(sizes); comp <- colnames(base$state); nC <- length(comp)
  d_time <- matrix(0, nS, length(traits), dimnames = list(NULL, traits))
  d_state <- array(0, c(nS, nC, length(traits)), dimnames = list(NULL, comp, traits))
  s0 <- individual$strategy$pars
  for (j in seq_along(traits)) {
    tr <- traits[j]
    d <- rel_step * max(abs(s0[[tr]]), 1e-8)
    rp <- perturb(tr, d); rm <- perturb(tr, -d)
    d_time[, j] <- (rp$time - rm$time) / (2 * d)
    d_state[, , j] <- (rp$state - rm$state) / (2 * d)
  }
  list(sizes = sizes, time = base$time, state = base$state,
       d_time = d_time, d_state = d_state)
}

# TF24f offspring_production trait gradient by reverse-mode AD (#472 scope B, the
# offspring surface). d(offspring_production)/d(theta) for the seed-rain integral
# offspring_production = sum_i tw_i * offspring_i (survival-weighted lifetime fecundity),
# the FROZEN rare-mutant / invasion gradient. The TF24 offspring tape with the one TF24f
# difference the census tape already solves -- the tracked collar is carried as a taped
# state with the curvature-linearised gradient-ascent rate (the lag the envelope theorem
# does not zero), seeded at the birth optimum (IFT-injected). Returns a named gradient
# vector + the reconstructed value.
tf24f_offspring_production_gradient <- function(scm, traits = NULL, species = 1L,
                                                birth_rate = NULL, trait_rel_step = 1e-5) {
  tf24f_require_strategy(scm)
  if (is.null(traits)) traits <- tf24_default_traits()
  # Fully native (mirrors FF16's offspring path): env + offspring harvest from the
  # live Patch; tf24f_harvest skipped. birth_rate < 0 recovers natively.
  strat <- scm$parameters$strategies[[species]]
  tf24f_offspring_gradient_native(scm, unlist(strat$pars), as.integer(species - 1L),
                                  if (is.null(birth_rate)) -1 else birth_rate,
                                  strat$k_acclim, strat$use_ad_gradient, traits,
                                  trait_rel_step)
}

# TF24f COUPLED census reconstruction -- the resident R0 gate (#472 scope B step 5,
# internal). A double-precision whole-stand re-evolution that reconstructs the ACTIVE
# canopy each RK stage and drives the real TF24f leaf at the reconstructed crown light.
# Returns the reconstructed TOTAL-stand census values and `env_err` (worst reconstructed
# vs SCM knot-light drift) -- the gauge that the coupled re-evolution reproduces the
# resident stand before the AD tape. Not a public gradient entry; used by the gate test.
tf24f_coupled_metrics <- function(scm, metrics = c("LAI", "size_moment"),
                                  species = 1L, birth_rate = NULL) {
  h <- tf24f_harvest(scm, species, birth_rate)
  if (length(h$nn_h) < 1L) {
    stop("the TF24f coupled (resident) census gradient needs the per-RK-stage ",
         "boundary-node harvest; re-run the resident SCM with ",
         "control(save_RK45_cache = TRUE)")
  }
  tf24f_coupled_metrics_impl(h$pp, h$eh, h$sh, h$birth_step, h$birth_rate, h$k_acclim,
                             h$use_ad_gradient, metrics, h$nn_h, h$nn_c, h$patch_area)
}

# TF24f RESIDENT (coupled) census trait gradient by reverse-mode AD (#472 scope B,
# build-order step 5 -- the resident-feedback AD refine of tf24f_resident_census_gradient_fd).
# The TOTAL stand-level d(census metric)/d(theta): all cohorts are re-evolved TOGETHER over
# the frozen schedule and the canopy light each RK stage is reconstructed from the active
# stand (heights AND densities respond to theta), so EVERY trait re-shades the light every
# cohort reads -- the genuine resident feedback that routinely flips the sign relative to
# the frozen (rare-mutant) gradient. The new ingredient beyond the frozen census tape is the
# leaf's crown-centre LIGHT channel (dprofit_dL + the collar-rate light curvature, harvested
# by FD over a flat environment); on the tape the resident light correction is anchored so
# the baseline is exact at theta0 and the focal plant's own height->light slope stays in the
# frozen-census dprofit_dh (no double count). One reverse sweep per metric over the coupled
# whole-stand replay. Validated against tf24f_resident_census_gradient_fd (full SCM re-runs);
# see test-tf24f-census-gradient.R. Single-species; FIXED node schedule (TF24f is stiff).
tf24f_resident_census_gradient_ad <- function(scm, metrics = c("LAI", "size_moment"),
                                              traits = NULL, species = 1L,
                                              birth_rate = NULL, trait_rel_step = 1e-5) {
  tf24f_require_strategy(scm)
  if (is.null(traits)) traits <- tf24_default_traits()
  # R0 gate (mirrors the multi-species path and FF16's coupled gate). The coupled double
  # re-evolution must reproduce the resident stand before the linearised sensitivity tape
  # is trusted. Past a short horizon the frozen-step replay drifts (joint env_err jumps
  # ~1e-6 -> 1e-2 beyond patch lifetime ~4) and the coupled log_density<->canopy
  # sensitivity runs away: TF24f's deeply-shaded leaf has a LARGE dprofit_dL, so the
  # high-density shaded cohorts strongly amplify the canopy-reshaping feedback -- a
  # stiffness the live SCM only tames by stepping ADAPTIVELY through it. The replay reuses
  # the SCM's frozen step sizes, so the linearised sensitivity diverges (finite-but-
  # astronomical, then NaN) with horizon. (FF16 stays robust to long horizons because its
  # closed-form net has a bounded light response -- no leaf-solve light amplification.)
  # Gate with a clear error rather than return NaN / garbage. feedback = "frozen" (the
  # invasion gradient) is robust at all horizons.
  diverged <- function(ee)
    stop("the TF24f coupled (resident) census re-evolution is too stiff on this node ",
         "schedule (joint env drift = ", signif(ee, 3), "). The single-species resident ",
         "census gradient is reliable only on a short, finely-resolved FIXED schedule ",
         "(patch lifetime ~4); use a shorter max_patch_lifetime or a finer fixed node ",
         "schedule, or use feedback = 'frozen' (the invasion gradient, robust at all ",
         "horizons).")
  r0 <- tf24f_coupled_metrics(scm, metrics = metrics, species = species,
                              birth_rate = birth_rate)
  if (!all(is.finite(r0$values)) || r0$env_err > 1e-2) diverged(r0$env_err)
  # Fully native gradient: env + birth steps + boundary-node history from the live Patch
  # (no tf24f_harvest, no Rcpp::as<> env). The boundary-node guard lives in the C++ entry.
  strat <- scm$parameters$strategies[[species]]
  g <- tf24f_coupled_gradient_native(scm, unlist(strat$pars), as.integer(species - 1L),
                                     if (is.null(birth_rate)) -1 else birth_rate,
                                     strat$k_acclim, strat$use_ad_gradient, traits, metrics,
                                     scm$parameters$patch_area, trait_rel_step)
  # Post-sweep guard: the linearised tape can still blow up (finite but astronomically
  # large) where the gate passes but the sensitivity is stiff; reject rather than mislead.
  if (!all(is.finite(g$jacobian)) || max(abs(g$jacobian)) > 1e12) diverged(r0$env_err)
  g
}

# Harvest a multi-species TF24f SCM into the all-species arrays the cross-species coupled
# engine consumes (#472 scope B, the cross-species resident Jacobian). The shared schedule
# + joint env, per-species parameter vectors / cohort birth steps / birth rates / TF24f
# acclimation knobs, and the all-species per-RK-stage boundary harvest
# (stand_newnode_*_stage_history_all, [step][stage][species]). The TF24f mirror of
# ff16_harvest_ms; the joint stand light is reconstructed in C++ from each species'
# re-evolved cohorts, so only these per-species pieces are needed from R.
tf24f_harvest_ms <- function(scm) {
  patch <- scm$patch
  nsp <- length(scm$parameters$strategies)
  sh  <- patch$step_history
  eh  <- patch$environment_history
  if (length(eh) < 1L) {
    stop("No resident schedule cached: run the SCM with control(save_RK45_cache = TRUE)")
  }
  nn_h <- patch$stand_newnode_height_stage_history_all
  nn_c <- patch$stand_newnode_competition_stage_history_all
  if (length(nn_h) < 1L) {
    stop("feedback = 'resident' on a multi-species TF24f stand needs the all-species ",
         "per-RK-stage harvest; re-run the resident SCM with control(save_RK45_cache = TRUE)")
  }
  list(
    pp_list = lapply(seq_len(nsp), function(s) unlist(scm$parameters$strategies[[s]]$pars)),
    eh = eh, sh = sh,
    birth_list = lapply(seq_len(nsp), function(s)
      vapply(patch$species[[s]]$node_times,
             function(t) which.min(abs(sh - t)) - 1L, integer(1))),
    birth_rate = vapply(seq_len(nsp), function(s)
      scm$offspring_production[[s]] / scm$net_reproduction_ratios[[s]], numeric(1)),
    k_acclim = vapply(seq_len(nsp), function(s)
      scm$parameters$strategies[[s]]$k_acclim, numeric(1)),
    use_ad_gradient = vapply(seq_len(nsp), function(s)
      as.integer(isTRUE(scm$parameters$strategies[[s]]$use_ad_gradient)), integer(1)),
    nn_h = nn_h, nn_c = nn_c,
    patch_area = scm$parameters$patch_area, nsp = nsp)
}

# TF24f CROSS-SPECIES resident (coupled) census trait gradient by reverse-mode AD (#472
# scope B). On a multi-species stand, the TOTAL-stand d(census metric)/d(theta of species
# `species`): all species' cohorts are re-evolved together over the frozen schedule and
# the canopy light each RK stage is the JOINT reconstruction (sum over species of each
# species' trapezium); only the target species' traits are registered as tape inputs, so
# one reverse sweep per metric gives the cross-species total -- the target's traits re-shade
# the joint canopy that EVERY species reads (the cross term the frozen gradient zeroes).
# A cheap double R0 pass gates a diverged (stiff) node schedule before the sweep. Unlike
# FF16 (closed-form leaf), the TF24f cross-species tape's linearised log_density / g'
# derivative is sensitive to the JOINT-canopy reshaping, and amplifies when the species
# are widely separated or the schedule is coarse -- so the gate is tighter (joint-env drift
# > 1e-3, vs FF16's 1e-2) and a post-sweep finiteness / magnitude guard catches any
# residual blow-up. The cross-species gradient therefore needs a WELL-CONDITIONED stand:
# closely-spaced species on a fine FIXED node schedule (where it matches the FD reference;
# see tf24f_resident_census_gradient_ms_fd). Validated on such a stand.
tf24f_resident_census_gradient_ms_ad <- function(scm, metrics = c("LAI", "size_moment"),
                                                 traits = NULL, species = 1L,
                                                 trait_rel_step = 1e-5) {
  tf24f_require_strategy(scm)
  if (is.null(traits)) traits <- tf24_default_traits()
  # Fully native (mirrors ff16_coupled_gradient_ms_native): the joint env + step schedule
  # + per-species birth steps + all-species boundary harvest come from the live Patch; only
  # the cheap per-species scalars (pp / recovered birth rates / acclimation knobs) are read
  # from $parameters in R. No tf24f_harvest_ms (no Rcpp::as<> env, no O(stand) patch rebuild).
  nsp <- length(scm$parameters$strategies)
  if (species < 1L || species > nsp) stop("target species out of range")
  pp_list <- lapply(seq_len(nsp), function(s) unlist(scm$parameters$strategies[[s]]$pars))
  br_vec <- vapply(seq_len(nsp), function(s)
    scm$offspring_production[[s]] / scm$net_reproduction_ratios[[s]], numeric(1))
  k_acclim <- vapply(seq_len(nsp), function(s)
    scm$parameters$strategies[[s]]$k_acclim, numeric(1))
  use_ad <- vapply(seq_len(nsp), function(s)
    as.integer(isTRUE(scm$parameters$strategies[[s]]$use_ad_gradient)), integer(1))
  patch_area <- scm$parameters$patch_area
  diverged <- function(ee)
    stop("the multi-species TF24f coupled re-evolution is too stiff on this node ",
         "schedule (joint env drift = ", signif(ee, 3), "). The cross-species gradient ",
         "needs a well-conditioned stand: closely-spaced species on a fine FIXED node ",
         "schedule (refine_schedule = FALSE). Re-run the resident SCM accordingly.")
  r0 <- tf24f_coupled_gradient_ms_native(scm, pp_list, br_vec, k_acclim, use_ad, traits,
          metrics, patch_area, as.integer(species), trait_rel_step, TRUE)
  if (!all(is.finite(r0$values)) || r0$env_err > 1e-3) diverged(r0$env_err)
  g <- tf24f_coupled_gradient_ms_native(scm, pp_list, br_vec, k_acclim, use_ad, traits,
         metrics, patch_area, as.integer(species), trait_rel_step, FALSE)
  # Post-sweep guard: the linearised tape can still blow up (finite but astronomically
  # large) where the gate passes but the derivative is stiff; reject rather than mislead.
  if (!all(is.finite(g$jacobian)) || max(abs(g$jacobian)) > 1e12) diverged(r0$env_err)
  g
}

# TF24f multi-species CROSS-SPECIES resident census FD reference (ground truth). Perturbs
# one species' (post-hyperpar) trait, re-runs the full SCM on the SAME fixed node schedule,
# and differences the TOTAL-stand metric (LAI = compute_competition(0); size_moment summed
# over species). Slow (one SCM solve per trait per side); keep `traits` small.
tf24f_resident_census_gradient_ms_fd <- function(p, env, ctrl,
                                                 metrics = c("LAI", "size_moment"),
                                                 traits = NULL, species = 1L,
                                                 rel_step = 1e-4) {
  if (is.null(traits)) traits <- tf24_default_traits()
  stand_metrics <- function(scm) {
    patch <- scm$patch
    out <- c(LAI = patch$compute_competition(0))
    sm <- 0
    for (s in seq_along(patch$species)) {
      sp <- patch$species[[s]]
      h <- sp$heights; dens <- exp(sp$log_densities)
      ord <- order(h, decreasing = TRUE); h <- h[ord]; dens <- dens[ord]
      phi <- dens * h
      if (length(h) > 1) sm <- sm + sum(0.5 * (h[-length(h)] - h[-1]) * (phi[-length(h)] + phi[-1]))
    }
    out["size_moment"] <- sm
    out[metrics]
  }
  run_with <- function(pp_override) {
    pmod <- p
    if (!is.null(pp_override)) {
      s <- pmod$strategies[[species]]; s$pars <- pp_override; pmod$strategies[[species]] <- s
    }
    run_scm(pmod, env, ctrl, refine_schedule = FALSE)
  }
  values <- stand_metrics(run_with(NULL))
  pars0 <- p$strategies[[species]]$pars
  jac <- matrix(0, length(metrics), length(traits), dimnames = list(metrics, traits))
  for (j in seq_along(traits)) {
    tr <- traits[j]
    if (!tr %in% names(pars0)) stop("unknown TF24f trait: ", tr)
    d <- rel_step * max(abs(pars0[[tr]]), 1e-8)
    pp_p <- pars0; pp_p[[tr]] <- pp_p[[tr]] + d
    pp_m <- pars0; pp_m[[tr]] <- pp_m[[tr]] - d
    jac[, j] <- (stand_metrics(run_with(pp_p)) - stand_metrics(run_with(pp_m))) / (2 * d)
  }
  list(jacobian = jac, values = values)
}

# TF24f RESIDENT (coupled) census trait gradient (#472 scope B, the resident-feedback
# surface -- prototype). The TOTAL stand-level d(census metric)/d(theta): unlike the
# frozen (rare-mutant) gradient, every cohort's height + density feeds back through the
# canopy light that the whole stand reads, so the feedback routinely dominates and can
# flip the sign relative to the frozen reading. Computed as a central finite difference
# over the FULL SCM: perturb a (post-hyperpar) strategy parameter, re-run run_scm on the
# SAME fixed node schedule, and difference the realised stand metric. This is the ground-
# truth resident gradient (the SCM responds in full) and the reference the coupled AD
# engine must reproduce; it is slow (one SCM solve per trait per side), so keep `traits`
# small. LAI is read from the realised patch (compute_competition(0)); size_moment is the
# size-distribution first moment Sum density_i * height_i (trapezium over the stand).
tf24f_resident_census_gradient_fd <- function(p, env, ctrl,
                                              metrics = c("LAI", "size_moment"),
                                              traits = NULL, species = 1L,
                                              rel_step = 1e-4) {
  if (is.null(traits)) traits <- tf24_default_traits()
  metric_set <- c("LAI", "size_moment")
  bad <- setdiff(metrics, metric_set)
  if (length(bad)) stop("unsupported resident metric(s): ", paste(bad, collapse = ", "))

  # Realised stand metrics from a completed SCM (the coupled, self-shaded stand).
  stand_metrics <- function(scm) {
    patch <- scm$patch
    out <- c(LAI = patch$compute_competition(0))
    sp <- patch$species[[species]]
    h <- sp$heights; dens <- exp(sp$log_densities)
    ord <- order(h, decreasing = TRUE)
    h <- h[ord]; dens <- dens[ord]
    # size_moment = trapezium of density*height over descending heights down to h0.
    phi <- dens * h
    sm <- 0
    if (length(h) > 1) sm <- sum(0.5 * (h[-length(h)] - h[-1]) * (phi[-length(h)] + phi[-1]))
    out["size_moment"] <- sm
    out[metrics]
  }
  run_with <- function(pp_override) {
    pmod <- p
    if (!is.null(pp_override)) {
      s <- pmod$strategies[[species]]
      s$pars <- pp_override
      pmod$strategies[[species]] <- s
    }
    run_scm(pmod, env, ctrl, refine_schedule = FALSE)
  }

  base_scm <- run_with(NULL)
  values <- stand_metrics(base_scm)
  pars0 <- p$strategies[[species]]$pars
  jac <- matrix(0, length(metrics), length(traits),
                dimnames = list(metrics, traits))
  for (j in seq_along(traits)) {
    tr <- traits[j]
    if (!tr %in% names(pars0)) stop("unknown TF24f trait: ", tr)
    d <- rel_step * max(abs(pars0[[tr]]), 1e-8)
    pp_p <- pars0; pp_p[[tr]] <- pp_p[[tr]] + d
    pp_m <- pars0; pp_m[[tr]] <- pp_m[[tr]] - d
    mp <- stand_metrics(run_with(pp_p)); mm <- stand_metrics(run_with(pp_m))
    jac[, j] <- (mp - mm) / (2 * d)
  }
  list(jacobian = jac, values = values)
}
