# Smoke test for the staging ATLS leaf demo (#566). Exercises the exact
# leaf-driving helpers the overstorey_staging/ vignette uses, so the demo cannot
# silently rot. overstorey_staging/ is .Rbuildignore'd (not installed), so this
# is dev-only: skip when the helper file is absent (installed/CRAN checks).

test_that("ATLS leaf demo helpers run end-to-end (damage, strategies, avoidance)", {
  skip_on_cran()
  helpers <- test_path("..", "..", "overstorey_staging", "atls_demo_helpers.R")
  skip_if_not(file.exists(helpers), "overstorey_staging/ not present (built package)")
  source(helpers, local = TRUE)

  # (1) The raw damage curve N(Tleaf): 1 when cold, in (0,1) when hot, monotone
  # non-increasing, finite everywhere.
  Ts <- seq(20, 50, by = 1)
  N <- atls_N_curve(Ts, atls_traits())
  expect_true(all(is.finite(N)))
  expect_true(all(N > 0 & N <= 1))
  expect_equal(N[Ts == 20], 1, tolerance = 1e-6)
  expect_true(N[Ts == 46] < 1)
  expect_true(all(diff(N) <= 1e-9))

  # (2) The four strategy axes move N in the predicted directions at a hot leaf.
  hot <- 42
  gen <- atls_N_curve(hot, atls_traits())
  expect_gt(atls_N_curve(hot, atls_traits(tcrit_0 = 44)), gen)   # tolerance
  expect_gt(atls_N_curve(hot, atls_traits(A_crit = 5)), gen)     # acclimation
  expect_gt(atls_N_curve(hot, atls_traits(k_r1_0 = 5000)), gen)  # repair

  # (3) Strategy solve gradient: every archetype x temperature solves finite,
  # with N in (0,1] and non-negative assimilation-side outputs sane.
  grad <- atls_solve_strategies(c(30, 38, 44), atls_strategies(), pm = FALSE)
  expect_true(all(vapply(grad$strategy, nzchar, logical(1))))
  for (v in c("Tleaf", "N", "A", "gs", "E", "profit")) {
    expect_true(all(is.finite(grad[[v]])), info = paste("non-finite", v))
  }
  expect_true(all(grad$N > 0 & grad$N <= 1))
  # A tolerant leaf keeps more carbon than the generalist at a hot leaf temp.
  tol44 <- grad$A[grad$strategy == "Tolerant"  & grad$Tenv == 44]
  gen44 <- grad$A[grad$strategy == "Generalist" & grad$Tenv == 44]
  expect_gt(tol44, gen44)

  # (4) Avoidance on the PM path: a better-coupled leaf runs cooler and retains
  # more capacity at the same air temperature.
  poor <- atls_solve_cell(1500, 34, 2, atls_traits(), pm = TRUE,
                          cfg = atls_leaf_config(d = 0.10,
                                                 leaf_specific_conductance_max = 2e-3))
  good <- atls_solve_cell(1500, 34, 2, atls_traits(), pm = TRUE,
                          cfg = atls_leaf_config(d = 0.02,
                                                 leaf_specific_conductance_max = 1e-2))
  expect_true(is.finite(poor$Tleaf) && is.finite(good$Tleaf))
  expect_lt(good$Tleaf, poor$Tleaf)   # cooler
  expect_gt(good$N, poor$N)           # -> higher N
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
