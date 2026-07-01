# Reverse-mode gradient of the emergent offspring_production for TF24 (#472 scope B,
# Phase F1-full). CI-runnable in plain R: the harvest + per-cohort reverse replay is
# compiled into plant.so (tf24_offspring_production_gradient_impl, the XAD adjoint tape
# resolved at load against odelia), so unlike the sourceCpp AD probes this needs no
# on-the-fly compilation.
#
# TF24 is STIFF (each rate eval re-solves the hydraulic leaf optimisation) and plants
# mature slowly (hmat ~ 16.6 m), so a faithful, fast CI stand is a tension: this uses a
# deliberately SMALL coarse-schedule stand (a handful of cohorts, short horizon) -- the
# reconstruction check is exact for ANY schedule, and lma's gradient (a growth/cascade
# trait) is robust even when offspring_production is tiny. Here we check the compiled API
# end-to-end: it reconstructs the SCM output and its gradient matches a two-pass FD for a
# representative trait (the per-trait FD sweep over the full trait set runs in the same
# style for any trait).

test_that("tf24_offspring_production_gradient reconstructs the SCM and matches FD", {
  H <- 6L                               # short horizon (stiff leaf opt is the cost)
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = 7))   # coarse -> few cohorts
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  scm <- run_scm(p, Environment("TF24"), ctlc, refine_schedule = FALSE)

  g <- offspring_production_gradient(scm)
  expect_length(g, 27L)                 # all net-production traits by default
  expect_true(all(is.finite(g)))
  expect_true("lma" %in% names(g) && "vcmax_25" %in% names(g))

  # (a) reconstruction: the harvest + reverse replay reproduce the SCM's emergent
  # output -- the strong end-to-end check (exercises the full kernel + tape).
  expect_equal(attr(g, "offspring_production"), scm$offspring_production[[1]],
               tolerance = 1e-3)

  # (b) gradient: two-pass FD over the SAME frozen schedule via the impl (re-harvests
  # the leaf opts + recomputes height_0 from the perturbed pp -> validates the cascade
  # + IFT-h0 terms). lma is the robust representative; tiny-magnitude photo traits are
  # vacuous at this small stand (validated at a real stand in the script).
  patch <- scm$patch          # cache: scm$patch rebuilds the whole patch each access
  sh <- patch$step_history
  eh <- patch$environment_history
  sp <- patch$species[[1]]
  nt <- sp$node_times
  pp <- unlist(scm$parameters$strategies[[1]]$pars)
  br <- scm$offspring_production[[1]] / scm$net_reproduction_ratios[[1]]
  birth_step <- vapply(nt, function(t) which.min(abs(sh - t)) - 1L, integer(1))
  N <- length(eh)
  tcoef <- numeric(length(nt)); x <- nt; n <- length(x)
  tcoef[1] <- 0.5 * (x[2] - x[1]); tcoef[n] <- 0.5 * (x[n] - x[n - 1])
  if (n > 2) tcoef[2:(n - 1)] <- 0.5 * (x[3:n] - x[1:(n - 2)])
  tw <- tcoef * sp$patch_densities * pp[["S_D"]] * br
  ah <- c(0, 0.2, 0.3, 0.6, 1.0, 0.875); hN <- diff(sh)
  ppsurv <- matrix(0, N, 6)
  for (k in seq_len(N)) for (s in 1:6) ppsurv[k, s] <- patch$pr_survival(sh[k] + ah[s] * hN[k])
  ppsab <- sp$pr_patch_survival_at_birth

  J_at <- function(v) {
    q <- pp; q[["lma"]] <- v
    attr(plant:::tf24_offspring_production_gradient_impl(q, eh, sh, birth_step, ppsurv,
                                                         ppsab, tw, "lma"),
         "offspring_production")
  }
  h <- 1e-4 * pp[["lma"]]
  fd_lma <- (J_at(pp[["lma"]] + h) - J_at(pp[["lma"]] - h)) / (2 * h)
  expect_equal(g[["lma"]], fd_lma, tolerance = 3e-3)
})

test_that("tf24_offspring_production_gradient selects the right species in a 2-sp stand", {
  # Two TF24 species sharing one canopy. The per-species harvest must select the
  # correct cohort family + strategy: each species' replay reconstructs ITS OWN
  # offspring_production[[s]]. (TF24 matures slowly, so at this small coarse stand the
  # emergent output is tiny -- reconstruction is exact for any magnitude; the gradient
  # is FD-validated in the single-species test above.)
  H <- 6L
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(c(0.1978791, 0.30), "lma"),
                      hyperpar = TF24_hyperpar, birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, H, length.out = 7), seq(0, H, length.out = 7))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  scm <- run_scm(p, Environment("TF24"), ctlc, refine_schedule = FALSE)
  expect_length(scm$offspring_production, 2L)

  for (s in 1:2) {
    gs <- offspring_production_gradient(scm, traits = "lma", species = s)
    expect_equal(attr(gs, "offspring_production"), scm$offspring_production[[s]],
                 tolerance = 1e-4)
  }
})
