## Unit tests for the scenario-evaluation engine (R/scenario_eval.R).
## These avoid long SCM runs except for one short smoke test.

test_that("scenario table and mapping load with expected structure", {
  tbl <- read_scenario_table()
  expect_s3_class(tbl, "tbl_df")
  expect_true(all(c("scenario_id", "Scenario", "Expectation",
                    "is_duplicate") %in% names(tbl)))
  expect_true(nrow(tbl) >= 1)
  ## The corrected CSV must be free of exact duplicate descriptor rows.
  expect_false(any(tbl$is_duplicate))

  map <- read_scenario_mapping()
  expect_true(all(c("csv_column", "level", "target", "param", "value") %in%
                    names(map)))
  expect_true(all(map$target %in% c("trait", "env", "driver")))
})

test_that("scenario_to_config translates a row into concrete settings", {
  tbl <- read_scenario_table()
  map <- read_scenario_mapping()
  cfg <- scenario_to_config(tbl[1, ], map)

  expect_true(cfg$expected %in% c("failure", "success"))
  ## K_s and lma are input traits; p_50/g1_TF24 are derived by TF24_hyperpar
  ## (#548) and must NOT be passed as input traits.
  expect_true(all(c("K_s", "lma", "rho", "vcmax_25", "theta",
                    "root_depth_shape_eta") %in% names(cfg$traits)))
  expect_false(any(c("p_50", "c", "b", "psi_crit", "g1_TF24") %in%
                     names(cfg$traits)))
  expect_true(is.numeric(cfg$traits) && all(is.finite(cfg$traits)))
  expect_true(!is.null(cfg$env$K_sat))
  expect_true(!is.null(cfg$driver$rainfall_mean))
})

test_that("scenario_to_config errors on an unmapped descriptor level", {
  tbl <- read_scenario_table()
  map <- read_scenario_mapping()
  bad <- tbl[1, ]
  bad$Ks <- "Nonsense level"
  expect_error(scenario_to_config(bad, map), "No mapping")
})

test_that("build_scenario applies traits and environment fields", {
  tbl <- read_scenario_table()
  map <- read_scenario_mapping()
  cfg <- scenario_to_config(tbl[1, ], map)
  built <- build_scenario(cfg, max_patch_lifetime = 5)

  pars <- built$p$strategies[[1]]$pars
  expect_equal(pars$lma, unname(cfg$traits[["lma"]]))
  expect_equal(pars$K_s, unname(cfg$traits[["K_s"]]))
  ## p_50 is derived from K_s by the hyperpar, so it must differ from the
  ## default once K_s is changed.
  expect_false(isTRUE(all.equal(pars$p_50, TF24_Strategy()$pars$p_50)))
  expect_equal(built$env$K_sat, cfg$env$K_sat)
})

test_that("classify_scm_run returns a well-formed classification", {
  ## A short, benign run using default traits and constant rainfall: should
  ## complete (not crash) regardless of demographic outcome.
  cfg <- list(traits = c(lma = 0.1978791), env = list(),
              driver = list(rainfall_mean = 1, rainfall_amp_frac = 0),
              expected = "success")
  built <- build_scenario(cfg, max_patch_lifetime = 2)
  run <- classify_scm_run(built$p, built$env, built$ctrl)

  expect_true(run$status %in% c("success", "failure"))
  expect_true(run$outcome %in% c("persisted", "extinct", "crashed"))
  expect_false(run$crashed)          # a benign run must not be a numerical crash
  expect_true(is.numeric(run$run_seconds))
})

test_that("scenario_model_fingerprint is a stable, non-empty hash", {
  fp <- scenario_model_fingerprint()
  expect_type(fp, "character")
  expect_true(nzchar(fp))
  expect_identical(fp, scenario_model_fingerprint())  # stable across calls
})

test_that("run_scenarios cache reuses unchanged scenarios", {
  ## Row 1 is mesic (LMA Low), row 3 is xeric (LMA High).
  scen <- read_scenario_table()[c(1, 3), ]
  map <- read_scenario_mapping()
  cache <- withr::local_tempfile(fileext = ".rds")

  cold <- run_scenarios(scen, map, max_patch_lifetime = 2, cache = cache)
  expect_true(file.exists(cache))

  ## A second identical run reuses everything and returns the same outcomes.
  expect_message(
    warm <- run_scenarios(scen, map, max_patch_lifetime = 2, cache = cache),
    "2 reused, 0 to run")
  expect_equal(warm$observed, cold$observed)

  ## Editing a mapping value used by one scenario invalidates only that one.
  map2 <- map
  map2$value[map2$csv_column == "LMA" & map2$level == "Low"] <- 0.05
  expect_message(
    run_scenarios(scen, map2, max_patch_lifetime = 2, cache = cache),
    "1 reused, 1 to run")

  ## A changed Control setting invalidates every scenario.
  expect_message(
    run_scenarios(scen, map, max_patch_lifetime = 2,
                  ctrl = control(ode_tol_rel = 1e-9), cache = cache),
    "0 reused, 2 to run")
})

test_that("scenario_summary tallies matches", {
  sc <- tibble::tibble(
    expected = c("failure", "success", "failure"),
    observed = c("failure", "success", "success"),
    match    = c(TRUE, TRUE, FALSE),
    persists = c(FALSE, TRUE, FALSE))
  s <- scenario_summary(sc)
  expect_equal(s$n, 3)
  expect_equal(s$n_match, 2)
  expect_equal(s$n_expected_fail_met, 1)
  expect_equal(s$n_expected_success_met, 1)
  ## Persistence is a separate axis from the match: a run can be a numerical
  ## success while the strategy dies out.
  expect_equal(s$n_persists, 1)
})

test_that("scenario_summary handles a scorecard recorded without persists", {
  ## The blessed baseline predates the column; summarising it must not error.
  sc <- tibble::tibble(
    expected = c("failure", "success"),
    observed = c("failure", "success"),
    match    = c(TRUE, TRUE))
  s <- scenario_summary(sc)
  expect_equal(s$n, 2)
  expect_equal(s$n_persists, 0)
})

test_that("persistence is judged at R0 >= 1, not R0 > 0", {
  ## The distinction the `persists` column exists to make: an offspring
  ## production of 2e-15 is extinction, not persistence, and the existing
  ## status/outcome pair calls it "persisted".
  expect_false(plant:::persists_at(2.2e-15, finite = TRUE))
  expect_false(plant:::persists_at(0.297, finite = TRUE))
  expect_true(plant:::persists_at(1, finite = TRUE))
  expect_true(plant:::persists_at(29.5, finite = TRUE))
  expect_false(plant:::persists_at(NA_real_, finite = FALSE))
})
