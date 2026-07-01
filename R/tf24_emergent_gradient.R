# TF24 emergent offspring_production gradient: there is no public tf24_*-named entry.
# offspring_production_gradient(scm) (R/emergent_gradient.R) dispatches on strategy and
# covers TF24 via stand_gradient()'s TF24 branch -> the native C++ entry
# tf24_offspring_production_gradient_native. (The single-strategy wrapper that used to
# live here was removed as redundant once stand_gradient unified the dispatch.)

# The 27 net-production TF24 trait names (10 leaf + 17 mass-cascade) the emergent
# gradients differentiate by default.
tf24_default_traits <- function() {
  c("vcmax_25","g1_TF24","beta2","K_s","b","c","jmax_25","a","curv_elec",
    "curv_colim","lma","rho","a_b1","r_l","r_b","r_s","r_r","k_l","k_b",
    "k_s","k_r","a_bio","a_y","a_l1","a_l2","theta","a_r1")
}

# Harvest a run-with-cache TF24 SCM into the frozen pieces the two-pass replay
# consumes (the ResidentHarvest seam, TF24 mirror of ff16_harvest). Requires the
# resident to have been run with shading_model = "crown-centre" (the replay re-solves
# the crown-centre leaf optimisation).
tf24_harvest <- function(scm, species = 1L, birth_rate = NULL) {
  types <- extract_RcppR6_template_types(scm$parameters, "Parameters")
  if (!identical(types[[1]], "TF24")) {
    stop("TF24 stand gradients are implemented for the TF24 strategy only")
  }
  if (species < 1L || species > length(scm$patch$species)) {
    stop("species index out of range: stand has ", length(scm$patch$species),
         " species")
  }
  # Cache the patch once: `scm$patch` rebuilds the whole RcppR6 patch object on each
  # access, so the per-stage pr_survival loop below is O(stand size) per call otherwise.
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
  pp    <- unlist(scm$parameters$strategies[[species]]$pars)

  if (is.null(birth_rate)) {
    birth_rate <- scm$offspring_production[[species]] / scm$net_reproduction_ratios[[species]]
  }

  # Cohort birth steps (introductions land on step times).
  birth_step <- vapply(nt, function(t) which.min(abs(sh - t)) - 1L, integer(1))
  N <- length(eh)
  # Node-spacing trapezoid weights so offspring_production == sum_i tw_i * offspring_i.
  tcoef <- numeric(length(nt)); x <- nt; n <- length(x)
  tcoef[1] <- 0.5 * (x[2] - x[1]); tcoef[n] <- 0.5 * (x[n] - x[n - 1])
  if (n > 2) tcoef[2:(n - 1)] <- 0.5 * (x[3:n] - x[1:(n - 2)])
  tw <- tcoef * pdens * pp[["S_D"]] * birth_rate
  # pr_patch_survival at the exact Cash-Karp stage times sh[k] + ah[s]*h.
  ah <- c(0, 0.2, 0.3, 0.6, 1.0, 0.875); hN <- diff(sh)
  ppsurv <- matrix(0, N, 6)
  for (k in seq_len(N)) for (s in 1:6) {
    ppsurv[k, s] <- patch$pr_survival(sh[k] + ah[s] * hN[k])
  }

  list(pp = pp, eh = eh, sh = sh, birth_step = birth_step, ppsurv = ppsurv,
       ppsab = ppsab, tw = tw, pdens = pdens, nt = nt, birth_rate = birth_rate)
}
