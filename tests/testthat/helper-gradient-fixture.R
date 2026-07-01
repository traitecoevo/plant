# Shared specification of the AD-gradient correctness-regression fixture
# (#472 scope B, refactor+optimize phase -- see notes/ad-refactor-optimize-roadmap.md).
#
# This is the "AD-vs-AD baseline": each spec computes a validated Jacobian on a
# fixed, deterministic stand and the fixture pins it to the committed snapshot
# (tests/testthat/fixtures/gradient-baseline.rds) so the upcoming refactor (moving
# the harvest into C++, unifying the three replay engines) is provably value-
# preserving. It is COMPLEMENTARY to the AD-vs-FD tests: those pin AD to the
# physics at ~1% (the recon noise floor); this pins AD to its own validated self
# at machine precision, catching a refactor regression in the 8th digit where FD
# cannot see.
#
# Tiers:
#   "bit"   -- pure-relocation / frozen-invasion paths -> bit-for-bit (1e-12 rel).
#   "noise" -- coupled / multi-species paths legitimately reorder FP sums (the
#              joint-canopy Sigma_species carries ~1e-6 summation noise) -> 5e-6 rel.
#
# Used by scripts/gradient_fixture.R (snapshot/check CLI) and
# tests/testthat/test-gradient-regression.R. Auto-sourced by testthat (helper-*.R)
# and source()d explicitly by the script after pkgload::load_all().

# ---- deterministic stands (fixed schedules; match the test fixtures) ----------

.gf_ff16_basic <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
.gf_ff16_resident <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 11)); p$max_patch_lifetime <- 60
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
.gf_ff16_ms <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, 70, length.out = 16), seq(0, 70, length.out = 16))
  p$max_patch_lifetime <- 70
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
.gf_ff16_statejac <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 9)); p$max_patch_lifetime <- 60
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
.gf_tf24f <- function(H = 4L, n = 9L) {
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = n))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
}
.gf_tf24 <- function(H = 8L, n = 9L) {
  p <- scm_base_parameters("TF24"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = n))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  run_scm(p, Environment("TF24"), ctlc, refine_schedule = FALSE)
}

.gf_tf24f_tr <- c("vcmax_25", "lma", "a_l1", "K_s")

# ---- the fixture specs ---------------------------------------------------------
# Each: compute() returns a list of named numeric arrays (jacobian/values/...);
# tier sets the comparison tolerance. Every distinct C++ engine the refactor
# touches appears exactly once.

gradient_fixture_specs <- function() {
  list(
    ff16_frozen_all = list(tier = "bit", compute = function() {
      scm <- .gf_ff16_basic()
      stand_gradient(scm, metrics = c("offspring_production", "LAI", "biomass",
                                      "size_moment"), traits = ff16_default_traits())
    }),
    ff16_offspring = list(tier = "bit", compute = function() {
      scm <- .gf_ff16_basic()
      list(gradient = offspring_production_gradient(scm, traits = ff16_default_traits()))
    }),
    ff16_resident = list(tier = "noise", compute = function() {
      scm <- .gf_ff16_resident()
      stand_gradient(scm, metrics = c("LAI", "biomass", "size_moment"),
                     traits = c("a_p1", "lma", "a_l1"), feedback = "resident")
    }),
    ff16_resident_ms = list(tier = "noise", compute = function() {
      scm <- .gf_ff16_ms()
      stand_gradient(scm, metrics = c("LAI", "size_moment"), traits = c("lma", "a_p1"),
                     species = 1L, feedback = "resident")
    }),
    ff16_birth_rate = list(tier = "noise", compute = function() {
      scm <- .gf_ff16_resident()
      g <- birth_rate_gradient(scm, metrics = c("LAI", "biomass", "size_moment",
                                                "net_reproduction_ratio"))
      list(d_birth_rate = g$d_birth_rate, values = g$values)
    }),
    ff16_state_jac = list(tier = "bit", compute = function() {
      scm <- .gf_ff16_statejac()
      J <- stand_state_jacobian(scm, traits = c("a_p1", "lma"))
      list(jacobian = J$jacobian, states = J$states)
    }),
    ff16_grow = list(tier = "bit", compute = function() {
      indv <- Individual("FF16", "FF16_Env")(FF16_Strategy())
      g <- grow_individual_to_size_gradient(indv, c(2, 5, 10), "height",
                                            Environment("FF16"), time_max = 200)
      list(d_state = g$d_state, d_time = g$d_time, time = as.numeric(g$time))
    }),
    tf24f_census = list(tier = "bit", compute = function() {
      scm <- .gf_tf24f()
      tf24f_census_gradient_ad(scm, metrics = c("LAI", "biomass", "size_moment"),
                               traits = .gf_tf24f_tr)
    }),
    tf24f_resident = list(tier = "noise", compute = function() {
      scm <- .gf_tf24f()
      tf24f_resident_census_gradient_ad(scm, metrics = c("LAI", "size_moment"),
                                        traits = .gf_tf24f_tr)
    }),
    tf24_offspring = list(tier = "bit", compute = function() {
      scm <- .gf_tf24()
      stand_gradient(scm, metrics = "offspring_production",
                     traits = c("vcmax_25", "lma", "a_l1", "theta"))
    })
  )
}

# Flatten a result list to a single named numeric vector for comparison.
gradient_fixture_flatten <- function(res) {
  keep <- intersect(c("jacobian", "values", "states", "d_state", "d_time", "time",
                      "gradient", "value"), names(res))
  out <- numeric(0)
  for (k in keep) {
    v <- res[[k]]
    if (is.list(v)) v <- unlist(v, use.names = FALSE)
    out <- c(out, stats::setNames(as.numeric(v),
                                  paste0(k, ".", seq_along(as.numeric(v)))))
  }
  out
}

gradient_fixture_tol <- function(tier) if (identical(tier, "bit")) 1e-12 else 5e-6
