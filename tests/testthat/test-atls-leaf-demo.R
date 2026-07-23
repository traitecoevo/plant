# Smoke test for the staging ATLS leaf demo (#566, merged revision). Exercises
# the exact leaf-driving helpers the overstorey_staging/ vignette uses, so the
# demo cannot silently rot. overstorey_staging/ is .Rbuildignore'd (not
# installed), so this is dev-only: skip when the helper file is absent.

test_that("ATLS leaf demo helpers run end-to-end (damage, strategies, avoidance)", {
  skip_on_cran()
  helpers <- test_path("..", "..", "overstorey_staging", "atls_demo_helpers.R")
  skip_if_not(file.exists(helpers), "overstorey_staging/ not present (built package)")
  source(helpers, local = TRUE)

  # (1) The deactivation substrate phi_d(Tleaf): ~0 when cold, rising through the
  # emergent onset (~34.5 C), bounded [0,1], monotone non-decreasing, finite.
  Ts <- seq(20, 50, by = 1)
  phi <- atls_phi_d_curve(Ts, atls_traits())
  expect_true(all(is.finite(phi)))
  expect_true(all(phi >= 0 & phi <= 1))
  expect_lt(phi[Ts == 20], 0.05)          # near-inert when cold
  expect_gt(phi[Ts == 45], 0.8)           # strongly unfolded when hot
  expect_true(all(diff(phi) >= -1e-9))    # monotone increasing in temperature

  # (2) The strategy axes lower standing recoverable damage I_r* at a hot leaf:
  # tolerance (shift onset up), acclimation (same, inducibly), repair (faster
  # resynthesis) each reduce I_r* below the generalist.
  hot <- 42
  gen <- atls_Ir_star_curve(hot, atls_traits())
  expect_lt(atls_Ir_star_curve(hot, atls_traits(topt_offset = 6)), gen)  # tolerance
  expect_lt(atls_Ir_star_curve(hot, atls_traits(A = 5)), gen)            # acclimation
  expect_lt(atls_Ir_star_curve(hot, atls_traits(k_rec = 2.0)), gen)      # repair

  # (3) Strategy solve gradient: every archetype x temperature solves finite, with
  # phi_d / I_r* bounded and assimilation-side outputs sane.
  grad <- atls_solve_strategies(c(30, 38, 44), atls_strategies(), pm = FALSE)
  expect_true(all(vapply(grad$strategy, nzchar, logical(1))))
  for (v in c("Tleaf", "phi_d", "Ir_star", "A", "gs", "E", "profit")) {
    expect_true(all(is.finite(grad[[v]])), info = paste("non-finite", v))
  }
  expect_true(all(grad$phi_d >= 0 & grad$phi_d <= 1))
  expect_true(all(grad$Ir_star >= 0 & grad$Ir_star <= 1))
  # A tolerant leaf keeps more carbon than the generalist at a hot leaf temp
  # (the thermostability offset raises the reversible capacity optimum).
  tol44 <- grad$A[grad$strategy == "Tolerant"   & grad$Tenv == 44]
  gen44 <- grad$A[grad$strategy == "Generalist" & grad$Tenv == 44]
  expect_gt(tol44, gen44)

  # (4) Avoidance on the PM path: a better-coupled leaf runs cooler and unfolds
  # less at the same air temperature.
  poor <- atls_solve_cell(1500, 34, 2, atls_traits(), pm = TRUE,
                          cfg = atls_leaf_config(d = 0.10,
                                                 leaf_specific_conductance_max = 2e-3))
  good <- atls_solve_cell(1500, 34, 2, atls_traits(), pm = TRUE,
                          cfg = atls_leaf_config(d = 0.02,
                                                 leaf_specific_conductance_max = 1e-2))
  expect_true(is.finite(poor$Tleaf) && is.finite(good$Tleaf))
  expect_lt(good$Tleaf, poor$Tleaf)   # cooler
  expect_lt(good$phi_d, poor$phi_d)   # -> less deactivation
})

test_that("ATLS demo SCM helpers run and rank strategies (community scale)", {
  skip_on_cran()
  helpers <- test_path("..", "..", "overstorey_staging", "atls_demo_helpers.R")
  skip_if_not(file.exists(helpers), "overstorey_staging/ not present (built package)")
  source(helpers, local = TRUE)

  # A single SCM run returns finite, sane community-fitness scalars.
  f <- atls_scm_fitness(25, type = "TF24t")
  expect_true(is.finite(f$R0) && f$R0 >= 0)
  expect_true(is.finite(f$offspring) && f$offspring >= 0)
  expect_identical(f$type, "TF24t")

  # Strategy tournament at the warm edge: tolerance persists best, and it beats
  # the generalist (the headline of the SCM section). One climate keeps it cheap.
  trn <- atls_scm_tournament(26)
  expect_setequal(trn$strategy,
                  c("Generalist", "Tolerant", "Repairer", "Acclimator"))
  expect_true(all(is.finite(trn$R0)))
  tol <- trn$R0[trn$strategy == "Tolerant"]
  gen <- trn$R0[trn$strategy == "Generalist"]
  expect_gt(tol, gen)
})
