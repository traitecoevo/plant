# AD-gradient correctness-regression fixture (#472 scope B, refactor+optimize phase).
# Pins each validated Jacobian to the committed snapshot
# (fixtures/gradient-baseline.rds) so the harvest->C++ move and the engine
# unification are provably value-preserving. Spec + setups: helper-gradient-fixture.R.
# Regenerate the snapshot after an INTENTIONAL change:
#   Rscript --no-init-file scripts/gradient_fixture.R snapshot
# Complementary to the AD-vs-FD tests (those pin AD to physics ~1%; this pins AD to
# its own validated self at machine precision).

test_that("AD gradients match the committed baseline (AD-vs-AD regression)", {
  # SAME-MACHINE guard. The baseline is snapshotted on one machine and the tiers are
  # bit (1e-12) / noise (5e-6) -- machine precision. Across compilers/libm the AD values
  # legitimately drift ~1e-9 (frozen, closed-form) to ~3e-3 (the coupled/resident paths,
  # whose cohort-height-crossing sort tie-breaks resolve differently), so committed FP
  # values cannot be asserted cross-platform. Skip on CI; the AD-vs-FD tests (loose,
  # physics-tolerance) are the portable correctness net. Run this on the snapshot machine
  # via `make test-ad` or `Rscript scripts/gradient_fixture.R check`.
  skip_on_ci()
  rds <- testthat::test_path("fixtures", "gradient-baseline.rds")
  skip_if_not(file.exists(rds),
              "no gradient-baseline.rds (run scripts/gradient_fixture.R snapshot)")
  # The AD engines are compiled into plant.so against the XAD tape; under a plain
  # load_all the FF16 tape symbols may be unresolved (see helper / odelia load order).
  skip_if_not(exists("ff16_stand_gradient_impl",
                     where = asNamespace("plant"), inherits = FALSE) ||
              is.function(tryCatch(plant:::ff16_stand_gradient_impl,
                                   error = function(e) NULL)),
              "AD tape symbols unavailable in this load_all session")

  baseline <- readRDS(rds)
  specs <- gradient_fixture_specs()

  for (nm in names(specs)) {
    tier <- specs[[nm]]$tier
    tol  <- gradient_fixture_tol(tier)
    ref  <- baseline[[nm]]
    expect_false(is.null(ref), info = paste0("no baseline entry for ", nm,
                 " -- regenerate the snapshot"))
    if (is.null(ref)) next
    cur <- gradient_fixture_flatten(specs[[nm]]$compute())
    expect_equal(length(cur), length(ref),
                 info = paste0(nm, ": fixture shape changed"))
    expect_equal(cur, ref, tolerance = tol,
                 info = paste0(nm, " (tier=", tier, ")"))
  }
})
