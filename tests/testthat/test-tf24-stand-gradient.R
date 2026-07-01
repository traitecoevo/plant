# Generic stand-gradient engine for TF24 (#472 scope B, build-order step 1, the TF24
# mirror). CI-runnable in plain R (the engine is compiled into plant.so, XAD tape
# resolved at load). TF24 is stiff (each rate re-solves the hydraulic leaf opt), so a
# deliberately small coarse-schedule stand keeps it fast; the gate + escape-hatch
# checks are exact for any schedule.

tf24_small_scm <- function(H = 8L, n = 9L) {
  p <- scm_base_parameters("TF24"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = n))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  run_scm(p, Environment("TF24"), ctlc, refine_schedule = FALSE)
}

test_that("stand_gradient (TF24) reproduces tf24_offspring_production_gradient", {
  scm <- tf24_small_scm()
  tr <- c("vcmax_25", "lma", "a_l1", "theta")
  g_ded <- as.numeric(offspring_production_gradient(scm, traits = tr))
  g_eng <- stand_gradient(scm, metrics = "offspring_production", traits = tr)
  expect_equal(as.numeric(g_eng$jacobian["offspring_production", ]), g_ded,
               tolerance = 1e-10)
  expect_equal(g_eng$values[["offspring_production"]], scm$offspring_production[[1]],
               tolerance = 1e-3)
})

test_that("stand_gradient (TF24) rejects census metrics (a documented follow-up)", {
  scm <- tf24_small_scm()
  expect_error(stand_gradient(scm, metrics = "LAI"), "follow-up")
  expect_error(stand_gradient(scm, metrics = "biomass"), "follow-up")
})

test_that("stand_state_jacobian (TF24) matches a frozen-resident FD per cohort", {
  scm <- tf24_small_scm()
  tr <- c("vcmax_25", "a_l1")
  J <- stand_state_jacobian(scm, traits = tr)
  nC <- nrow(J$states)
  expect_equal(dim(J$jacobian), c(nC, 6L, length(tr)))

  # Frozen-resident FD of the cohort final states, perturbing vcmax_25 (a leaf trait
  # that does not move height_0, so the FD is not limited by the seedling root-find).
  h <- plant:::tf24_harvest(scm, 1L, NULL)
  states_at <- function(q) plant:::tf24_state_jacobian_impl(
    q, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, tr, h$birth_rate)$states
  hh <- 1e-6 * abs(h$pp[["vcmax_25"]]); q1 <- q2 <- h$pp
  q1[["vcmax_25"]] <- q1[["vcmax_25"]] + hh; q2[["vcmax_25"]] <- q2[["vcmax_25"]] - hh
  fd <- (states_at(q1) - states_at(q2)) / (2 * hh)
  ti <- match("vcmax_25", tr)
  for (i in c(2L, nC - 1L)) for (cc in c("height", "fecundity")) {
    ci <- match(cc, dimnames(J$states)[[2]])
    expect_equal(unname(J$jacobian[i, ci, ti]), unname(fd[i, ci]), tolerance = 1e-4)
  }
})
