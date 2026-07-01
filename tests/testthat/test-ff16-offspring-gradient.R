# Reverse-mode gradient of the emergent offspring_production (#472 scope B, Phase C).
# CI-runnable in plain R: the reverse-mode replay is compiled into plant.so
# (ff16_offspring_production_gradient_impl, the XAD adjoint tape resolved at load
# against odelia), so unlike the sourceCpp AD probes this needs no on-the-fly
# compilation and runs everywhere.

test_that("offspring_production_gradient matches a two-pass finite difference", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)

  g <- offspring_production_gradient(scm, traits = c("a_p1", "lma"))
  # lma changes the seedling size height_0 (a height_seed root-find), so its gradient
  # exercises the implicit-function-theorem h0 path, not just the demographic replay.

  # The replay reconstructs the SCM's emergent output.
  expect_equal(attr(g, "offspring_production"), scm$offspring_production[[1]],
               tolerance = 1e-4)

  # Two-pass FD over the SAME frozen schedule: gather the replay inputs, perturb the
  # trait in the parameter vector, finite-difference the reconstructed value.
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

  # Two-pass FD over the same frozen schedule, perturbing one trait in the parameter
  # vector. The impl recomputes height_0 from the (perturbed) parameters, so this is
  # an h0-ACTIVE finite difference -- it validates the IFT seedling-size term too.
  fd_best <- function(trait, rel_h) {
    J_at <- function(v) {
      q <- pp; q[[trait]] <- v
      attr(ff16_offspring_production_gradient_impl(q, eh, sh, birth_step, ppsurv,
                                                   ppsab, tw, trait),
           "offspring_production")
    }
    fds <- vapply(rel_h, function(rh) {
      h <- rh * abs(pp[[trait]])
      (J_at(pp[[trait]] + h) - J_at(pp[[trait]] - h)) / (2 * h)
    }, numeric(1))
    fds[which.min(abs(fds - g[[trait]]))]   # best step (FD has a truncation/roundoff sweet spot)
  }
  # a_p1: physiology, does not touch height_0.
  expect_equal(g[["a_p1"]], fd_best("a_p1", c(1e-5, 1e-6)), tolerance = 1e-4)
  # lma: flows through the IFT height_0 term as well as the replay.
  expect_equal(g[["lma"]], fd_best("lma", c(1e-4, 1e-5)), tolerance = 1e-3)
})

test_that("offspring_production_gradient handles a multi-species stand per species", {
  # Two FF16 species (distinct lma) sharing one canopy. offspring_production is a
  # per-species vector; the resident light (environment_history) is the SHARED frozen
  # canopy of BOTH species, so each species' cohorts replay against the same env --
  # the gradient is the rare-mutant / invasion gradient of species `s` against the
  # fixed two-species canopy. Coarse schedule keeps it fast.
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2178), "lma"),
                      hyperpar = FF16_hyperpar, birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 9), seq(0, 60, length.out = 9))
  p$max_patch_lifetime <- 60
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  expect_length(scm$offspring_production, 2L)

  # (a) Each species reconstructs ITS OWN offspring_production[[s]] -- the strong
  # end-to-end check that the per-species harvest selects the right cohort family.
  for (s in 1:2) {
    gs <- offspring_production_gradient(scm, traits = "a_p1", species = s)
    expect_equal(attr(gs, "offspring_production"), scm$offspring_production[[s]],
                 tolerance = 1e-4)
  }

  # (b) The species-2 gradient matches a two-pass FD over species 2's frozen schedule.
  g2 <- offspring_production_gradient(scm, traits = "a_p1", species = 2L)
  patch <- scm$patch          # cache: scm$patch rebuilds the whole patch each access
  sh <- patch$step_history
  eh <- patch$environment_history
  sp <- patch$species[[2]]
  nt <- sp$node_times
  pp <- unlist(scm$parameters$strategies[[2]]$pars)
  br <- scm$offspring_production[[2]] / scm$net_reproduction_ratios[[2]]
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
    q <- pp; q[["a_p1"]] <- v
    attr(plant:::ff16_offspring_production_gradient_impl(q, eh, sh, birth_step, ppsurv,
                                                         ppsab, tw, "a_p1"),
         "offspring_production")
  }
  h <- 1e-6 * pp[["a_p1"]]
  fd <- (J_at(pp[["a_p1"]] + h) - J_at(pp[["a_p1"]] - h)) / (2 * h)
  expect_equal(g2[["a_p1"]], fd, tolerance = 1e-4)
})
