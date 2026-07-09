## Opt-in gateway check: does the current build still reproduce the recorded
## per-scenario outcomes? This is a *baseline diff*, not an "all pass" assertion
## — many scenarios are expected to fail by design. A changed outcome (a
## regression, or an improvement once NSC/#548 land) fails the test so it can be
## consciously re-blessed by regenerating the recorded baseline with
## `make bless-scenarios` (writes tests/testthat/test_data/scenario_baseline.rds).
##
## Slow (runs the full SCM for every scenario); gated behind an env var.

test_that("scenario outcomes match the recorded baseline", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("PLANT_RUN_SCENARIOS")),
              "Set PLANT_RUN_SCENARIOS=1 to run the (slow) scenario gateway.")

  baseline_path <- test_path("test_data", "scenario_baseline.rds")
  skip_if_not(file.exists(baseline_path), "No recorded scenario baseline.")
  baseline <- readRDS(baseline_path)

  mpl <- attr(baseline, "max_patch_lifetime")
  if (is.null(mpl)) mpl <- 100
  workers <- as.integer(Sys.getenv("SCENARIO_WORKERS", "1"))
  current <- run_scenarios(max_patch_lifetime = mpl, workers = workers)

  b <- baseline[order(baseline$scenario_id), ]
  cur <- current[order(current$scenario_id), ]
  expect_equal(cur$scenario_id, b$scenario_id)

  ## Compare the binary success/failure classification only: offspring
  ## magnitudes can differ across platforms, but the classification is stable.
  changed <- b$observed != cur$observed
  if (any(changed)) {
    msg <- paste(sprintf("%s: %s -> %s", b$scenario_id[changed],
                         b$observed[changed], cur$observed[changed]),
                 collapse = "; ")
    fail(paste0("Scenario outcomes changed vs baseline ",
                "(re-bless via scripts/run_scenario_gateway.R if intended): ", msg))
  } else {
    succeed()
  }
})
