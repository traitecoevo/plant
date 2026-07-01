# TF24f stand census reconstruction -- the R0 GATE (#472 scope B, build-order step 1).
# TF24f is the right target for the stand CENSUS
# gradients (its tracked-collar leaf evaluation is analytic / IFT-able, so the census
# number density needs no curvature harvest). Before any reverse-mode tape (R1) is built,
# this gate proves the prerequisite §7 flags: that re-evolving each cohort's {5 demog,
# tracked collar, log_density} over the SCM's frozen Cash-Karp schedule -- driving the
# real TF24f leaf at the tracked collar in double precision -- reproduces the SCM's stored
# stand (heights / opt_root_psi_state / log_density) and the census LAI =
# compute_competition(0). The recon is value-only (no AD), so it runs in plain R / CI.
#
# Patch lifetime 4 keeps the tallest cohort ~7 m, inside the regime where the replay is
# bit-faithful (heights / collar to ~1e-6, log_density to ~1e-6). For taller cohorts the
# growth-rate-gradient finite difference (g' = backward FD of height_dt, abs step 1e-6) is
# the fidelity-limiting term: as growth saturates with height, g' is a near-cancelling
# difference of two large crown-centre leaf solves, and at lifetime 5 (tallest ~8.4 m) the
# oldest cohort's log_density drifts ~6e-3 (LAI ~0.4%); at lifetime 8 the tracked collar
# drifts out of the leaf transport spline domain entirely (an acclimation
# schedule-sensitivity to pin before pushing to long horizons).

tf24f_small_scm <- function(H = 4L, n = 9L) {
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = n))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
}

# Read the SCM's stored per-cohort final states. node_at() is 1-based; the ODE state
# layout is [height, mortality, fecundity, area_heartwood, mass_heartwood,
# opt_root_psi_state, offspring, log_density].
tf24f_scm_stand <- function(scm, species = 1L) {
  sp <- scm$patch$species[[species]]
  nc <- length(sp$node_times)
  list(
    heights = vapply(seq_len(nc), function(i) sp$node_at(i)$ode_state[1], numeric(1)),
    collar  = vapply(seq_len(nc), function(i) sp$node_at(i)$ode_state[6], numeric(1)),
    log_density = sp$log_densities,
    LAI = scm$patch$compute_competition(0))
}

test_that("tf24f_census_recon reproduces the SCM's collar-state stand (R0 gate)", {
  scm <- tf24f_small_scm()
  ref <- tf24f_scm_stand(scm)
  rec <- plant:::tf24f_census_recon(scm, metrics = c("LAI", "biomass", "size_moment"))

  # Per-cohort heights and the tracked collar (opt_root_psi_state) -- the replay
  # re-evolves the 6-state ODE, so matching these confirms the collar-state trajectory.
  expect_equal(rec$heights, ref$heights, tolerance = 1e-4)
  expect_equal(rec$collar,  ref$collar,  tolerance = 1e-4)
  # The census number density (the genuinely new state). g' (the SCM's backward-FD
  # growth-rate gradient) is the fidelity-limiting term, hence the looser absolute floor.
  expect_equal(rec$log_density, ref$log_density, tolerance = 1e-3)
})

test_that("tf24f_census_recon LAI reduction matches compute_competition(0)", {
  scm <- tf24f_small_scm()
  ref <- tf24f_scm_stand(scm)
  rec <- plant:::tf24f_census_recon(scm, metrics = "LAI")
  # The height-trapezium of density * k_I * area_leaf (+ pending-seed tail) is exactly
  # the SCM's competition integral at ground level.
  expect_equal(rec$values[["LAI"]], ref$LAI, tolerance = 1e-3)
  # biomass / size_moment are the same trapezium with a different per-cohort kernel;
  # they must be finite and positive on a live stand.
  rec2 <- plant:::tf24f_census_recon(scm, metrics = c("biomass", "size_moment"))
  expect_true(all(is.finite(rec2$values)) && all(rec2$values > 0))
})

test_that("tf24f_census_gradient_fd: finite frozen census Jacobian (R1 gate)", {
  # The reverse-mode census trait gradient (R1) must reproduce a two-pass FD over the
  # SAME frozen-env replay; this builds that FD reference and checks it is well-formed
  # and self-consistent (a usable prototype gradient pending the AD tape).
  scm <- tf24f_small_scm()
  tr <- c("vcmax_25", "lma", "a_l1", "K_s")
  g <- plant:::tf24f_census_gradient_fd(scm, metrics = c("LAI", "biomass", "size_moment"),
                                        traits = tr)
  expect_equal(dim(g$jacobian), c(3L, length(tr)))
  expect_true(all(is.finite(g$jacobian)))
  # The reconstructed metric values match the recon (and so the SCM census).
  expect_equal(g$values[["LAI"]], scm$patch$compute_competition(0), tolerance = 1e-3)
  # Self-consistency: a column recomputed at a 10x coarser FD step agrees to ~1% -- the
  # precision of any FD over the recon, capped by the leaf optimiser/root-find noise
  # floor (~5e-7) amplified by the step. This is the tolerance band the reverse-mode AD
  # tape will be validated within (the AD removes this FD-step sensitivity).
  g2 <- plant:::tf24f_census_gradient_fd(scm, metrics = "LAI", traits = "lma",
                                         rel_step = 1e-4)
  expect_equal(unname(g$jacobian["LAI", "lma"]), unname(g2$jacobian["LAI", "lma"]),
               tolerance = 2e-2)
})

test_that("tf24f_census_gradient_ad reproduces the frozen census FD gradient (R1 tape)", {
  # The reverse-mode AD census tape (build-order step 2) is the REFINE of the per-trait
  # FD over the recon: ONE reverse sweep over the 7-state replay {5 demog, tracked collar,
  # log_density} per metric. The collar is carried as a taped state with a curvature-
  # linearised gradient-ascent rate (the envelope theorem does NOT zero its theta-
  # contribution -- it lags the optimum at k_acclim = 1), and the census g' reproduces the
  # SCM's backward-FD scheme via a second harvested operating point at h - GEPS. It must
  # match tf24f_census_gradient_fd to ~1% (the recon noise floor); the AD removes the FD
  # step-size sensitivity, so it is the more accurate of the two.
  scm <- tf24f_small_scm()
  tr <- c("vcmax_25", "lma", "a_l1", "K_s")
  mt <- c("LAI", "biomass", "size_moment")
  ad <- plant:::tf24f_census_gradient_ad(scm, metrics = mt, traits = tr)
  fd <- plant:::tf24f_census_gradient_fd(scm, metrics = mt, traits = tr)

  expect_equal(dim(ad$jacobian), c(length(mt), length(tr)))
  expect_true(all(is.finite(ad$jacobian)))
  # The reconstructed metric VALUES are the recon's (the linearisation is exact at theta0),
  # so they match the SCM census and the FD gate bit-for-bit.
  expect_equal(unname(ad$values[["LAI"]]), scm$patch$compute_competition(0), tolerance = 1e-3)
  expect_equal(unname(ad$values), unname(fd$values), tolerance = 1e-8)
  # Gradient agreement: AD vs FD within the recon noise floor (~1%); use a relative
  # comparison so the comparison is not dominated by the large-magnitude entries.
  rel <- abs(ad$jacobian - fd$jacobian) / pmax(abs(fd$jacobian), 1e-8)
  expect_lt(max(rel), 0.02)
})

test_that("stand_gradient dispatches TF24f census to the AD tape", {
  # The first-class API (build-order step 3): stand_gradient(strat = TF24f) routes the
  # census metrics under feedback = "frozen" to the reverse-mode AD tape, and refuses the
  # not-yet-built surfaces (resident coupling, offspring_production) with a clear message.
  scm <- tf24f_small_scm()
  tr <- c("vcmax_25", "lma", "a_l1", "K_s")
  mt <- c("LAI", "biomass", "size_moment")
  sg <- stand_gradient(scm, metrics = mt, traits = tr)
  ad <- plant:::tf24f_census_gradient_ad(scm, metrics = mt, traits = tr)
  expect_equal(sg$jacobian, ad$jacobian)
  expect_equal(sg$values, ad$values)
  expect_equal(unname(sg$values[["LAI"]]), scm$patch$compute_competition(0), tolerance = 1e-3)
  # feedback = "resident" routes to the coupled total gradient (step 5).
  sgr <- stand_gradient(scm, metrics = "LAI", traits = tr, feedback = "resident")
  adr <- plant:::tf24f_resident_census_gradient_ad(scm, metrics = "LAI", traits = tr)
  expect_equal(sgr$jacobian, adr$jacobian)
  # offspring_production routes to the offspring tape, and assembles alongside census.
  sgo <- stand_gradient(scm, metrics = c("offspring_production", "LAI"), traits = tr)
  go <- plant:::tf24f_offspring_production_gradient(scm, traits = tr)
  expect_equal(unname(sgo$jacobian["offspring_production", ]), as.numeric(go$gradient))
  expect_equal(unname(sgo$jacobian["LAI", ]), unname(ad$jacobian["LAI", ]))
})

test_that("tf24f offspring_production AD == FD over the same reconstruction", {
  # The offspring tape (the TF24 offspring tape + the tracked-collar curvature the census
  # tape already solves) must exactly differentiate its frozen-invasion reconstruction.
  # The rigorous check is AD vs a central FD over the SAME offspring reconstruction
  # (re-harvested at perturbed traits), as for the census tape; and the reconstructed value
  # matches the SCM's offspring_production.
  scm <- tf24f_small_scm()
  tr <- c("vcmax_25", "lma", "a_l1", "K_s")
  ad <- plant:::tf24f_offspring_production_gradient(scm, traits = tr)
  expect_equal(ad$value, scm$offspring_production[[1]], tolerance = 1e-3)
  h <- plant:::tf24f_harvest(scm)
  oval <- function(pp) plant:::tf24f_offspring_gradient_impl(pp, h$eh, h$sh, h$birth_step,
    h$ppsurv, h$ppsab, h$tw, h$k_acclim, h$use_ad_gradient, tr, 1e-5)$value
  fd <- vapply(tr, function(t) {
    d <- 1e-5 * max(abs(h$pp[[t]]), 1e-8)
    pp_p <- h$pp; pp_p[[t]] <- pp_p[[t]] + d
    pp_m <- h$pp; pp_m[[t]] <- pp_m[[t]] - d
    (oval(pp_p) - oval(pp_m)) / (2 * d)
  }, numeric(1))
  rel <- abs(ad$gradient - fd) / pmax(abs(fd), 1e-30)
  expect_lt(max(rel), 0.02)   # the recon noise floor (most entries agree to ~1e-5)
})

# A well-conditioned 2-species TF24f stand (closely-spaced lma, fine fixed schedule) -- the
# regime where the cross-species coupled tape is stable (the joint-canopy reshaping does not
# blow up the linearised log_density / g' derivative).
tf24f_ms_scm <- function(lmas = c(0.15, 0.22), n = 12L, H = 4L) {
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(lmas, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = as.list(rep(20, length(lmas))))
  p$node_schedule_times <- rep(list(seq(0, H, length.out = n)), length(lmas))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  list(p = p, ctlc = ctlc, scm = run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE))
}

test_that("tf24f multi-species coupled R0 reconstructs the joint canopy + total metrics", {
  w <- tf24f_ms_scm()
  hm <- plant:::tf24f_harvest_ms(w$scm)
  tr <- c("lma", "K_s"); mt <- c("LAI", "size_moment")
  r0 <- plant:::tf24f_coupled_gradient_ms_impl(hm$pp_list, hm$eh, hm$sh, hm$birth_list,
    hm$birth_rate, hm$k_acclim, hm$use_ad_gradient, tr, mt, hm$nn_h, hm$nn_c,
    hm$patch_area, 1L, 1e-5, TRUE)
  expect_lt(r0$env_err, 1e-3)
  expect_equal(unname(r0$values[["LAI"]]), w$scm$patch$compute_competition(0), tolerance = 1e-3)
})

test_that("tf24f cross-species coupled AD == FD over the same recon; cross term + gate", {
  # The cross-species total gradient d(total-stand metric)/d(theta of species 1): all
  # species' cohorts re-evolved together, the joint canopy active in species 1's traits.
  # Validated AD == FD over the SAME coupled reconstruction on LAI/lma (the dominant entry,
  # as for the FF16 ms test), and the resident cross-species gradient differs by orders of
  # magnitude from the frozen rare-mutant gradient (the canopy feedback).
  w <- tf24f_ms_scm()
  hm <- plant:::tf24f_harvest_ms(w$scm)
  tr <- c("lma", "K_s")
  ad <- stand_gradient(w$scm, metrics = "LAI", traits = tr, species = 1L, feedback = "resident")
  recon <- function(ppl) plant:::tf24f_coupled_gradient_ms_impl(ppl, hm$eh, hm$sh,
    hm$birth_list, hm$birth_rate, hm$k_acclim, hm$use_ad_gradient, tr, "LAI", hm$nn_h,
    hm$nn_c, hm$patch_area, 1L, 1e-5, TRUE)$values
  d <- 1e-5 * abs(hm$pp_list[[1]][["lma"]])
  ppp <- hm$pp_list; ppp[[1]][["lma"]] <- ppp[[1]][["lma"]] + d
  ppm <- hm$pp_list; ppm[[1]][["lma"]] <- ppm[[1]][["lma"]] - d
  fd_lma <- (recon(ppp) - recon(ppm)) / (2 * d)
  expect_equal(unname(ad$jacobian["LAI", "lma"]), unname(fd_lma[["LAI"]]), tolerance = 0.05)
  expect_equal(unname(ad$values[["LAI"]]), w$scm$patch$compute_competition(0), tolerance = 1e-3)
  # The cross-species feedback dominates: resident differs hugely from frozen.
  gf <- plant:::tf24f_census_gradient_ad(w$scm, metrics = "LAI", traits = "lma", species = 1L)
  expect_false(isTRUE(all.equal(unname(gf$jacobian["LAI", "lma"]),
                                unname(ad$jacobian["LAI", "lma"]), tolerance = 0.5)))
  # An ill-conditioned stand (widely-spaced species, coarse schedule) is rejected, not
  # silently returned as a blown-up gradient.
  b <- tf24f_ms_scm(lmas = c(0.1218, 0.2625), n = 8L)
  expect_error(stand_gradient(b$scm, metrics = "LAI", traits = tr, species = 1L,
                              feedback = "resident"), "well-conditioned")
})

test_that("tf24f coupled census R0 reconstructs the resident stand (step 5 gate)", {
  # The resident (coupled) gate: a double-precision whole-stand re-evolution that
  # reconstructs the ACTIVE canopy each RK stage (every cohort re-shades the light the
  # whole stand reads) must reproduce the SCM's resident stand before the AD tape. The
  # gauge is env_err (worst reconstructed vs SCM knot-light drift) + the LAI value.
  scm <- tf24f_small_scm()
  r0 <- plant:::tf24f_coupled_metrics(scm, metrics = c("LAI", "size_moment"))
  expect_lt(r0$env_err, 1e-4)
  expect_equal(unname(r0$values[["LAI"]]), scm$patch$compute_competition(0), tolerance = 1e-3)
  expect_true(all(is.finite(r0$values)) && all(r0$values > 0))
})

test_that("tf24f coupled census establishes the final-boundary cohort (a_l2 NaN fix)", {
  # Regression (mirrors the FF16 fix): the last node's introduction time lands on the
  # final step (birth step == N), and the coupled re-evolution loop (rn < N) never
  # established it -- so its census height stayed zero-initialised, biasing the value
  # and making the a_l2 channel of tf24_area_leaf hit 0*log(0) = NaN (which poisons the
  # whole shared census tape). The boundary cohort must now be established at h0.
  scm <- tf24f_small_scm()
  # The coupled VALUE reconstructs the frozen-engine value (both include the boundary
  # cohort at h0) and stays a sane positive census.
  cm <- plant:::tf24f_coupled_metrics(scm, metrics = c("LAI", "biomass", "size_moment"))
  expect_true(all(is.finite(cm$values)) && all(cm$values > 0))
  # The resident gradient is finite for a trait set INCLUDING a_l2 (the allometry
  # exponent that triggered the 0*log(0) NaN through area_leaf = (h/a_l1)^(1/a_l2)).
  tr <- c("vcmax_25", "lma", "a_l1", "a_l2", "K_s", "theta")
  g <- plant:::tf24f_resident_census_gradient_ad(
    scm, metrics = c("LAI", "biomass", "size_moment"), traits = tr)
  expect_true(all(is.finite(g$jacobian)))
  expect_true(is.finite(g$jacobian["LAI", "a_l2"]))
})

test_that("tf24f resident coupled AD == FD over the same coupled reconstruction (step 5)", {
  # The coupled AD tape (build-order step 5) must exactly differentiate the coupled
  # whole-stand reconstruction it replays. Mirroring the FF16 coupled test, the rigorous
  # check is AD vs a central FD over the SAME reconstruction (tf24f_coupled_metrics_impl)
  # -- isolating tape correctness from the looser full-SCM gaps (the adaptive node
  # schedule's trait response and metric-definition differences). The new ingredient over
  # the frozen tape -- the leaf crown-centre LIGHT channel + anchored resident reshaping --
  # is what makes EVERY trait feed back through the canopy.
  scm <- tf24f_small_scm()
  tr <- c("lma", "K_s")
  mt <- c("LAI", "size_moment")
  h <- plant:::tf24f_harvest(scm)
  recon <- function(pp) plant:::tf24f_coupled_metrics_impl(pp, h$eh, h$sh, h$birth_step,
    h$birth_rate, h$k_acclim, h$use_ad_gradient, mt, h$nn_h, h$nn_c, h$patch_area)$values
  ad <- plant:::tf24f_resident_census_gradient_ad(scm, metrics = mt, traits = tr)
  fd <- matrix(0, length(mt), length(tr), dimnames = list(mt, tr))
  for (j in seq_along(tr)) {
    d <- 1e-5 * max(abs(h$pp[[tr[j]]]), 1e-8)
    pp_p <- h$pp; pp_p[[tr[j]]] <- pp_p[[tr[j]]] + d
    pp_m <- h$pp; pp_m[[tr[j]]] <- pp_m[[tr[j]]] - d
    fd[, j] <- (recon(pp_p) - recon(pp_m)) / (2 * d)
  }
  # AD values match the coupled reconstruction to its floor (~env_err): the AD tape
  # linearises the leaf at the frozen env, while tf24f_coupled_metrics drives the real
  # leaf at the RECONSTRUCTED canopy light -- the two differ only by the knot-light drift.
  expect_equal(unname(ad$values), unname(recon(h$pp)), tolerance = 1e-5)
  # AD == FD over the coupled reconstruction within the recon noise floor.
  rel <- abs(ad$jacobian - fd) / pmax(abs(fd), 1e-8)
  expect_lt(max(rel), 0.02)
})

test_that("tf24f resident LAI gradient flips sign and tracks the full-SCM FD", {
  # The defining resident property: the canopy feedback flips d(LAI)/d(lma) negative
  # (a stand-wide LMA rise shades out leaf area) where the frozen rare-mutant gradient is
  # positive. The coupled AD reproduces this and lands near the full-SCM FD ground truth
  # for LAI (the canonical compute_competition(0) metric); a residual gap is expected from
  # the adaptive node schedule's own trait response (the frozen-schedule replay holds it
  # fixed), as for the FF16 coupled gradient.
  H <- 4L
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = 9L))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  scm <- run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)

  gf <- plant:::tf24f_census_gradient_ad(scm, metrics = "LAI", traits = "lma")
  gr <- plant:::tf24f_resident_census_gradient_ad(scm, metrics = "LAI", traits = "lma")
  fd <- plant:::tf24f_resident_census_gradient_fd(p, Environment("TF24f"), ctlc,
          metrics = "LAI", traits = "lma")
  expect_gt(gf$jacobian["LAI", "lma"], 0)         # frozen: positive (rare mutant)
  expect_lt(gr$jacobian["LAI", "lma"], 0)         # resident: feedback flips it negative
  expect_lt(fd$jacobian["LAI", "lma"], 0)         # full-SCM FD agrees on the sign
  # Same order of magnitude as the full-SCM ground truth (grid-response gap aside).
  expect_equal(unname(gr$jacobian["LAI", "lma"]),
               unname(fd$jacobian["LAI", "lma"]), tolerance = 0.2)
})

test_that("tf24f resident census FD gradient is the total gradient (flips sign)", {
  # The resident (coupled) gradient is the TOTAL stand-level d(metric)/d(theta): the
  # canopy feedback routinely dominates and flips the sign relative to the frozen
  # (rare-mutant / invasion) gradient. Ground truth here is FD over the full SCM.
  H <- 4L
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = 9L))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  tr <- c("lma", "K_s")
  gr <- plant:::tf24f_resident_census_gradient_fd(p, Environment("TF24f"), ctlc,
          metrics = c("LAI", "size_moment"), traits = tr)
  expect_equal(dim(gr$jacobian), c(2L, length(tr)))
  expect_true(all(is.finite(gr$jacobian)))

  scm <- run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
  expect_equal(unname(gr$values[["LAI"]]), scm$patch$compute_competition(0),
               tolerance = 1e-6)
  # The defining property: feedback flips d(LAI)/d(lma) -- positive for a rare mutant
  # against the frozen canopy, negative as the whole stand's LMA rises.
  gf <- plant:::tf24f_census_gradient_fd(scm, metrics = "LAI", traits = "lma")
  expect_gt(gf$jacobian["LAI", "lma"], 0)
  expect_lt(gr$jacobian["LAI", "lma"], 0)
})

test_that("tf24f census recon is strategy-guarded", {
  # A TF24 (non-f) resident has no tracked-collar state; the harvest must reject it.
  p <- scm_base_parameters("TF24"); p$max_patch_lifetime <- 4L
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 4L, length.out = 9L))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  scm <- run_scm(p, Environment("TF24"), ctlc, refine_schedule = FALSE)
  expect_error(plant:::tf24f_census_recon(scm), "TF24f strategy only")
})

test_that("tf24f single-species coupled gradient gates the long-horizon stiff replay", {
  # The coupled (resident) census gradient re-evolves the whole stand over the SCM's
  # FROZEN step schedule. Past a short horizon (~patch lifetime 4) that frozen-step
  # replay drifts and the coupled log_density<->canopy sensitivity runs away -- TF24f's
  # deeply-shaded leaf has a large dprofit_dL, so high-density shaded cohorts amplify a
  # canopy-reshaping feedback the live SCM only tames by stepping ADAPTIVELY. The
  # linearised sensitivity then diverges (finite-but-astronomical -> NaN) with horizon.
  # The path must GATE that with a clear error (mirroring the multi-species path and
  # FF16's coupled gate), not return NaN / garbage.
  scm4 <- tf24f_small_scm(H = 4L)
  g <- plant:::tf24f_resident_census_gradient_ad(scm4, metrics = c("LAI", "size_moment"),
                                                 traits = c("vcmax_25", "lma"))
  expect_true(all(is.finite(g$jacobian)))          # short horizon: still reliable
  # A long horizon is rejected with an actionable message, NOT a NaN/astronomical Jacobian.
  scm6 <- tf24f_small_scm(H = 6L)
  expect_error(
    plant:::tf24f_resident_census_gradient_ad(scm6, metrics = c("LAI", "size_moment"),
                                              traits = c("vcmax_25", "lma")),
    "too stiff")
  # The public entry gates identically.
  expect_error(stand_gradient(scm6, metrics = "LAI", traits = "lma", feedback = "resident"),
               "too stiff")
  # The FROZEN (invasion) gradient stays robust at the same long horizon.
  gf <- plant:::tf24f_census_gradient_ad(scm6, metrics = c("LAI", "size_moment"),
                                         traits = c("vcmax_25", "lma"))
  expect_true(all(is.finite(gf$jacobian)))
})
