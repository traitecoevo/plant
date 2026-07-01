# Generic calibration-facing stand-gradient engine (#472 scope B, build-order
# step 1). CI-runnable in plain R: the engine (ff16_stand_gradient_impl,
# ff16_state_jacobian_impl) is compiled into plant.so with the XAD adjoint tape
# resolved at load against odelia, so no on-the-fly compilation is needed.

test_that("stand_gradient reproduces offspring_production_gradient as one entry", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)

  tr <- c("a_p1", "lma", "theta", "omega")
  g_ded <- as.numeric(offspring_production_gradient(scm, traits = tr))
  g_eng <- stand_gradient(scm, metrics = "offspring_production", traits = tr)
  # The generic engine is bit-for-bit the dedicated routine for this metric.
  expect_equal(as.numeric(g_eng$jacobian["offspring_production", ]), g_ded,
               tolerance = 1e-10)
  expect_equal(g_eng$values[["offspring_production"]], scm$offspring_production[[1]],
               tolerance = 1e-4)
})

test_that("stand_gradient census metrics reconstruct + match a frozen-resident FD", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)

  mets <- c("offspring_production", "LAI", "biomass", "size_moment")
  tr   <- c("a_p1", "lma", "a_l1")
  g <- stand_gradient(scm, metrics = mets, traits = tr)

  # LAI reconstructs the SCM's own compute_competition(0) (the size-distribution
  # leaf-area integral, including the pending-seed tail term).
  expect_equal(g$values[["LAI"]], scm$patch$compute_competition(0), tolerance = 1e-6)
  expect_equal(g$values[["offspring_production"]], scm$offspring_production[[1]],
               tolerance = 1e-4)
  expect_equal(dim(g$jacobian), c(length(mets), length(tr)))

  # Frozen-resident two-pass FD of each metric value (perturb a trait in the
  # parameter vector, re-reduce over the SAME harvested schedule + resident light).
  h <- plant:::ff16_harvest(scm, 1L, NULL)
  val_at <- function(q) plant:::ff16_stand_gradient_impl(
    q, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, tr, mets, h$birth_rate,
    "frozen", list(), list(), h$patch_area, -1, -1)$values
  fd_best <- function(trait, mm, rel_h) {
    fds <- vapply(rel_h, function(rh) {
      hh <- rh * abs(h$pp[[trait]]); q1 <- q2 <- h$pp
      q1[[trait]] <- q1[[trait]] + hh; q2[[trait]] <- q2[[trait]] - hh
      (val_at(q1)[[mm]] - val_at(q2)[[mm]]) / (2 * hh)
    }, numeric(1))
    fds[which.min(abs(fds - g$jacobian[mm, trait]))]
  }
  for (mm in mets) {
    expect_equal(g$jacobian[mm, "a_p1"], fd_best("a_p1", mm, c(1e-5, 1e-6)),
                 tolerance = 1e-3)
    # lma exercises the IFT seedling-size (height_0) path through every metric.
    expect_equal(g$jacobian[mm, "lma"], fd_best("lma", mm, c(1e-4, 1e-5)),
                 tolerance = 1e-3)
  }
})

test_that("feedback='resident' is the coupled total gradient (every trait feeds back)", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 11))
  p$max_patch_lifetime <- 60
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)

  mets <- c("LAI", "biomass", "size_moment")
  tr   <- c("a_p1", "lma", "a_l1")
  gr <- stand_gradient(scm, metrics = mets, traits = tr, feedback = "resident")
  gf <- stand_gradient(scm, metrics = mets, traits = tr, feedback = "frozen")
  expect_equal(dim(gr$jacobian), c(length(mets), length(tr)))

  # The coupled feedback makes a_p1 and lma (NOT just the allometric a_l1) move the
  # canopy -> the resident gradient differs materially from the frozen one.
  expect_true(all(abs(gr$jacobian - gf$jacobian) > 0))

  # AD vs a central FD over the SAME coupled reconstruction (frozen geometry). The
  # reconstruction carries ~1e-8 value noise, so validate at the noise-optimal step.
  h <- plant:::ff16_harvest(scm, 1L, NULL)
  cval <- function(q) plant:::ff16_coupled_metrics_impl(
    q, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, mets, h$birth_rate,
    h$nn_h, h$nn_c, h$patch_area)$values
  for (trait in c("a_p1", "lma")) {
    d <- 1e-4 * abs(h$pp[[trait]]); q1 <- q2 <- h$pp
    q1[[trait]] <- q1[[trait]] + d; q2[[trait]] <- q2[[trait]] - d
    fd <- (cval(q1) - cval(q2)) / (2 * d)
    for (mm in mets)
      expect_equal(unname(gr$jacobian[mm, trait]), unname(fd[[mm]]),
                   tolerance = 5e-3)
  }

  # offspring_production stays the FROZEN invasion gradient even under "resident".
  gro <- stand_gradient(scm, metrics = "offspring_production", traits = tr,
                        feedback = "resident")
  gfo <- stand_gradient(scm, metrics = "offspring_production", traits = tr,
                        feedback = "frozen")
  expect_equal(gro$jacobian, gfo$jacobian, tolerance = 1e-8)
})

test_that("birth_rate_gradient is the resident-coupled d(census)/d(birth_rate)", {
  # A new gradient AXIS (the birth-rate driver, not a trait): birth_rate scales every
  # cohort's establishment density, re-shading the whole resident canopy. The
  # resident-coupled derivative is the demographic-equilibrium object; it is genuinely
  # non-trivial (the canopy feedback dominates and even flips signs), unlike the frozen
  # reading which is the trivial identity d(metric)/d(birth_rate) = metric/birth_rate.
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 11))
  p$max_patch_lifetime <- 60
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  mets <- c("LAI", "biomass", "size_moment")

  g <- birth_rate_gradient(scm, metrics = mets)
  expect_equal(names(g$d_birth_rate), mets)
  expect_true(all(is.finite(g$d_birth_rate)))
  # LAI reconstruction matches the SCM's own compute_competition(0) (recon faithfulness).
  expect_equal(unname(g$values[["LAI"]]), scm$patch$compute_competition(0),
               tolerance = 1e-3)

  # AD == central FD over the SAME coupled reconstruction (perturb the birth_rate driver
  # of ff16_coupled_metrics_impl, everything else held). The coupled recon carries ~1e-8
  # value noise, so validate at the coupled noise floor (~1e-2), as the trait test does.
  h <- plant:::ff16_harvest(scm, 1L, NULL)
  cval_br <- function(br) plant:::ff16_coupled_metrics_impl(
    h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, mets, br,
    h$nn_h, h$nn_c, h$patch_area)$values
  br0 <- h$birth_rate; d <- 1e-5 * br0
  fd <- (cval_br(br0 + d) - cval_br(br0 - d)) / (2 * d)
  for (mm in mets)
    expect_equal(unname(g$d_birth_rate[[mm]]), unname(fd[[mm]]), tolerance = 1e-2)

  # The resident feedback makes it genuinely different from the trivial frozen identity
  # (metric / birth_rate) -- the whole point of differentiating the coupled stand. LAI
  # differs by >50% relative; biomass even flips sign (the self-shading feedback dominates).
  frozen_identity <- g$values / br0
  expect_gt(abs(g$d_birth_rate[["LAI"]] - frozen_identity[["LAI"]]) /
              abs(frozen_identity[["LAI"]]), 0.5)
  expect_true(g$d_birth_rate[["biomass"]] * frozen_identity[["biomass"]] < 0)

  # offspring_production's birth-rate derivative is the trivial frozen identity -> excluded.
  expect_error(birth_rate_gradient(scm, metrics = "offspring_production"),
               "census metrics")
})

test_that("birth_rate_gradient net_reproduction_ratio = dR0/db (the R0=1 / mutant axis)", {
  # The demographic-equilibrium axis: d(net_reproduction_ratio)/d(birth_rate). With unit
  # per-seed weights, birth_rate enters R0 ONLY through the resident canopy (denser
  # recruitment -> more self-shading -> lower per-seed reproduction), so this is the pure
  # density feedback -- equivalently, the change in a trait-identical mutant's fitness as
  # the resident density moves. It is what an R0 = 1 self-replacement Newton step needs.
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 11))
  p$max_patch_lifetime <- 60
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  g <- birth_rate_gradient(scm, metrics = "net_reproduction_ratio")
  # The reconstructed R0 reproduces the SCM's own net_reproduction_ratio.
  expect_equal(unname(g$values[["net_reproduction_ratio"]]),
               scm$net_reproduction_ratios[[1]], tolerance = 1e-6)
  # Density-dependent: dR0/db < 0 (the feedback that makes an R0 = 1 equilibrium well-posed).
  expect_lt(g$d_birth_rate[["net_reproduction_ratio"]], 0)

  # AD == central FD over the SAME coupled recon, holding the per-seed weight fixed at the
  # harvest birth_rate (birth_rate0) and perturbing only the canopy/density driver -- the
  # mutant reading. (Perturbing the weight too would add a spurious 1/b term.)
  h <- plant:::ff16_harvest(scm, 1L, NULL)
  cval <- function(br) plant:::ff16_coupled_metrics_impl(h$pp, h$eh, h$sh, h$birth_step,
    h$ppsurv, h$ppsab, h$tw, "net_reproduction_ratio", br, h$nn_h, h$nn_c, h$patch_area,
    TRUE, h$birth_rate)$values
  br0 <- h$birth_rate; d <- 1e-5 * br0
  fd <- (cval(br0 + d) - cval(br0 - d)) / (2 * d)
  expect_equal(unname(g$d_birth_rate[["net_reproduction_ratio"]]), unname(fd[[1]]),
               tolerance = 1e-2)

  # net_reproduction_ratio is the single-species demographic-equilibrium axis (cross-
  # species is the census metrics only).
  pms <- scm_base_parameters("FF16")
  pms <- add_strategies(pms, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                        birth_rate = list(20, 20))
  pms$node_schedule_times <- list(seq(0, 60, length.out = 11), seq(0, 60, length.out = 11))
  pms$max_patch_lifetime <- 60
  scm_ms <- run_scm(pms, Environment("FF16"), control(save_RK45_cache = TRUE),
                    refine_schedule = FALSE)
  expect_error(birth_rate_gradient(scm_ms, metrics = "net_reproduction_ratio"),
               "single-species")
})

# Helper: harvest a multi-species FF16 SCM into the per-species arrays the multi-species
# coupled engine consumes (mirrors ff16_harvest's single-species gather, all species).
ms_harvest <- function(scm) {
  patch <- scm$patch
  nsp <- length(scm$parameters$strategies)
  sh  <- patch$step_history
  list(
    pp_list = lapply(seq_len(nsp), function(s)
      unlist(scm$parameters$strategies[[s]]$pars)),
    eh = patch$environment_history, sh = sh,
    birth_list = lapply(seq_len(nsp), function(s)
      vapply(patch$species[[s]]$node_times,
             function(t) which.min(abs(sh - t)) - 1L, integer(1))),
    birth_rate = vapply(seq_len(nsp), function(s)
      scm$offspring_production[[s]] / scm$net_reproduction_ratios[[s]], numeric(1)),
    nn_h = patch$stand_newnode_height_stage_history_all,
    nn_c = patch$stand_newnode_competition_stage_history_all,
    area = scm$parameters$patch_area, nsp = nsp)
}

test_that("feedback='resident' on a multi-species stand is the cross-species gradient", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, 70, length.out = 16), seq(0, 70, length.out = 16))
  p$max_patch_lifetime <- 70
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  tr <- c("lma", "a_p1")
  # The public resident path on a >1-species stand routes to the cross-species coupled
  # engine: d(TOTAL-stand metric)/d(theta of species `species`), finite + matching the
  # internal impl, and materially different from the per-species frozen reading.
  gr <- stand_gradient(scm, metrics = c("LAI", "size_moment"), traits = tr,
                       species = 1L, feedback = "resident")
  expect_true(all(is.finite(gr$jacobian)))
  h <- plant:::ff16_harvest_ms(scm)
  gms <- plant:::ff16_coupled_gradient_ms_impl(h$pp_list, h$eh, h$sh, h$birth_list, tr,
           c("LAI", "size_moment"), h$birth_rate, h$nn_h, h$nn_c, h$patch_area, 1L)
  expect_equal(gr$jacobian, gms$jacobian, tolerance = 1e-10)
  gf <- stand_gradient(scm, metrics = c("LAI", "size_moment"), traits = tr,
                       species = 1L, feedback = "frozen")
  expect_true(all(abs(gr$jacobian - gf$jacobian) > 0))
  # frozen is also fine multi-species (per-species invasion gradient).
  expect_equal(dim(gf$jacobian), c(2L, length(tr)))
})

test_that("multi-species coupled R0 reconstructs the joint canopy + total metrics", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, 70, length.out = 14), seq(0, 70, length.out = 14))
  p$max_patch_lifetime <- 70
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  h <- ms_harvest(scm)
  mets <- c("LAI", "biomass", "size_moment")
  ms <- plant:::ff16_coupled_metrics_ms_impl(h$pp_list, h$eh, h$sh, h$birth_list, mets,
          h$birth_rate, h$nn_h, h$nn_c, h$area)
  # joint canopy re-evolution tracks the SCM stand (env drift small on a fixed schedule)
  expect_lt(ms$env_err, 1e-3)
  # total LAI reconstructs the SCM's compute_competition(0) at the final census
  scm_lai <- -log(h$eh[[length(h$eh)]][[6]]$get_environment_at_height(0))
  expect_equal(unname(ms$values[["LAI"]]), scm_lai, tolerance = 5e-3)
  # total biomass/size_moment == the single-species frozen engine summed over species
  frz <- c(LAI = 0, biomass = 0, size_moment = 0)
  for (s in seq_len(h$nsp))
    frz <- frz + stand_gradient(scm, metrics = mets, species = s,
                                feedback = "frozen")$values
  expect_equal(unname(ms$values[["biomass"]]), unname(frz[["biomass"]]),
               tolerance = 5e-3)
})

test_that("multi-species cross-species gradient: AD == coupled-recon FD, cross term != 0", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, 70, length.out = 16), seq(0, 70, length.out = 16))
  p$max_patch_lifetime <- 70
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  h <- ms_harvest(scm)
  mets <- c("LAI", "size_moment"); tr <- "a_p1"; target <- 1L
  g <- plant:::ff16_coupled_gradient_ms_impl(h$pp_list, h$eh, h$sh, h$birth_list, tr,
         mets, h$birth_rate, h$nn_h, h$nn_c, h$area, target)
  # AD vs central FD over the SAME coupled reconstruction. Validate on LAI, the
  # cleanest metric: the mass/size-weighted reductions carry cohort-height-crossing
  # value noise that swamps a coarse-schedule FD (a step sweep shows size_moment's FD
  # bounces by 10x while AD is stable -- AD is the exact reference, FD the noise floor).
  cval <- function(pt) {
    pl <- h$pp_list; pl[[target]] <- pt
    plant:::ff16_coupled_metrics_ms_impl(pl, h$eh, h$sh, h$birth_list, mets,
      h$birth_rate, h$nn_h, h$nn_c, h$area)$values
  }
  d <- 1e-4 * abs(h$pp_list[[target]][[tr]])
  p1 <- p2 <- h$pp_list[[target]]; p1[[tr]] <- p1[[tr]] + d; p2[[tr]] <- p2[[tr]] - d
  fd <- (cval(p1) - cval(p2)) / (2 * d)
  expect_equal(unname(g$jacobian["LAI", tr]), unname(fd[["LAI"]]), tolerance = 2e-2)
  # the cross-species + self feedback moves the total gradient off the frozen reading
  gf <- stand_gradient(scm, metrics = mets, traits = tr, species = target,
                       feedback = "frozen")
  expect_true(all(abs(g$jacobian[, tr] - gf$jacobian[, tr]) > 0))
})

test_that("birth_rate_gradient cross-species: AD == coupled-recon FD (species 1 driver)", {
  # The cross-species birth-rate axis: d(TOTAL-stand metric)/d(birth_rate of species 1) on
  # a 2-species fixed-schedule stand. Species 1's recruitment density re-shades the joint
  # canopy every species reads -- the genuinely new cross term (the frozen reading is the
  # diagonal identity with zero cross term). Validated AD == FD over the SAME coupled recon.
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, 70, length.out = 16), seq(0, 70, length.out = 16))
  p$max_patch_lifetime <- 70
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  mets <- c("LAI", "size_moment")
  g <- birth_rate_gradient(scm, metrics = mets, species = 1L)
  expect_true(all(is.finite(g$d_birth_rate)))
  expect_equal(unname(g$values[["LAI"]]), scm$patch$compute_competition(0), tolerance = 1e-3)

  h <- ms_harvest(scm)
  cval <- function(brv) plant:::ff16_coupled_metrics_ms_impl(h$pp_list, h$eh, h$sh,
    h$birth_list, mets, brv, h$nn_h, h$nn_c, h$area)$values
  br0 <- h$birth_rate; d <- 1e-5 * br0[1]
  brp <- brm <- br0; brp[1] <- brp[1] + d; brm[1] <- brm[1] - d
  fd <- (cval(brp) - cval(brm)) / (2 * d)
  expect_equal(unname(g$d_birth_rate[["LAI"]]), unname(fd[["LAI"]]), tolerance = 2e-2)
})

test_that("coupled resident gradient is finite for the FULL default trait set (a_l2 NaN fix)", {
  # Regression: a cohort whose introduction time lands on the final step (birth step == N)
  # was never established in the coupled re-evolution, so its census height stayed at the
  # zero-initialised value and the a_l2 channel of area_leaf = (h/a_l1)^(1/a_l2) hit
  # 0*log(0) = NaN (and the corrupted zero-height cohort also biased the census VALUE).
  # Single- and multi-species, all 28 default traits must now be finite.
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 70, length.out = 16)); p$max_patch_lifetime <- 70
  scm1 <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                  refine_schedule = FALSE)
  g1 <- stand_gradient(scm1, metrics = c("LAI", "biomass", "size_moment"),
                       traits = ff16_default_traits(), feedback = "resident")
  expect_true(all(is.finite(g1$jacobian)))
  expect_true(is.finite(g1$jacobian["LAI", "a_l2"]) && g1$jacobian["LAI", "a_l2"] != 0)

  p2 <- scm_base_parameters("FF16")
  p2 <- add_strategies(p2, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                       birth_rate = list(20, 20))
  p2$node_schedule_times <- list(seq(0, 70, length.out = 16), seq(0, 70, length.out = 16))
  p2$max_patch_lifetime <- 70
  scm2 <- run_scm(p2, Environment("FF16"), control(save_RK45_cache = TRUE),
                  refine_schedule = FALSE)
  g2 <- stand_gradient(scm2, metrics = c("LAI", "size_moment"),
                       traits = ff16_default_traits(), species = 1L, feedback = "resident")
  expect_true(all(is.finite(g2$jacobian)))
})

test_that("multi-species (> 2 species) coupled resident: finite, all targets, total == frozen-summed", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.15, 0.25), "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20, 20, 20))
  p$node_schedule_times <- rep(list(seq(0, 70, length.out = 14)), 3)
  p$max_patch_lifetime <- 70
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  expect_equal(length(scm$parameters$strategies), 3L)
  mets <- c("LAI", "biomass", "size_moment")
  # Every target species' cross-species total gradient is finite for all 28 traits.
  for (tgt in 1:3) {
    g <- stand_gradient(scm, metrics = mets, traits = ff16_default_traits(),
                        species = tgt, feedback = "resident")
    expect_true(all(is.finite(g$jacobian)))
  }
  # R0: the coupled TOTAL-stand value reconstructs the frozen-engine sum over species.
  h <- ms_harvest(scm)
  ms <- plant:::ff16_coupled_metrics_ms_impl(h$pp_list, h$eh, h$sh, h$birth_list, mets,
          h$birth_rate, h$nn_h, h$nn_c, h$area)
  frz <- c(LAI = 0, biomass = 0, size_moment = 0)
  for (s in 1:3) frz <- frz + stand_gradient(scm, metrics = mets, species = s,
                                             feedback = "frozen")$values
  expect_equal(unname(ms$values), unname(frz[mets]), tolerance = 5e-3)
})

test_that("stand_state_jacobian matches a frozen-resident FD per cohort", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 9))
  p$max_patch_lifetime <- 60
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)

  tr <- c("a_p1", "lma")
  J <- stand_state_jacobian(scm, traits = tr)
  nC <- nrow(J$states)
  expect_equal(dim(J$jacobian), c(nC, 6L, length(tr)))

  # Frozen-resident FD of the cohort final states, perturbing a_p1.
  h <- plant:::ff16_harvest(scm, 1L, NULL)
  states_at <- function(q) plant:::ff16_state_jacobian_impl(
    q, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, tr)$states
  hh <- 1e-6 * abs(h$pp[["a_p1"]]); q1 <- q2 <- h$pp
  q1[["a_p1"]] <- q1[["a_p1"]] + hh; q2[["a_p1"]] <- q2[["a_p1"]] - hh
  fd <- (states_at(q1) - states_at(q2)) / (2 * hh)
  ti <- match("a_p1", tr)
  # Check a mid and a late cohort's height + offspring sensitivity.
  for (i in c(2L, nC - 1L)) for (cc in c("height", "offspring")) {
    ci <- match(cc, dimnames(J$states)[[2]])
    expect_equal(unname(J$jacobian[i, ci, ti]), unname(fd[i, ci]),
                 tolerance = 1e-3)
  }
})
