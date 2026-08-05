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

  ## Compare the two *classifications*, not the magnitudes: offspring production
  ## differs across platforms, but which side of a threshold it falls on is
  ## stable.
  ##
  ## `observed` alone is not enough. It tests `finite && total > 0`, which every
  ## scenario has satisfied since the crash fixes landed, so it is pinned at
  ## 8/8 "success" and cannot move. Measured: switching the density coordinate
  ## (#590) multiplies R0 by 2.4x to 47x and moves S01 across R0 = 1, and
  ## `observed` does not change on a single scenario. A guard that survives
  ## that is not guarding the thing the gateway now reports, so `persists` --
  ## the R0 >= 1 axis #572 made the headline -- is diffed as well.
  report <- function(col, label) {
    if (!col %in% names(b) || !col %in% names(cur)) {
      return(NULL)
    }
    changed <- b[[col]] != cur[[col]]
    if (!any(changed, na.rm = TRUE)) {
      return(NULL)
    }
    sprintf("%s %s", label,
            paste(sprintf("%s: %s -> %s", b$scenario_id[changed],
                          b[[col]][changed], cur[[col]][changed]),
                  collapse = "; "))
  }

  msgs <- c(report("observed", "observed"), report("persists", "persists"))
  if (length(msgs)) {
    fail(paste0("Scenario outcomes changed vs baseline ",
                "(re-bless via `make bless-scenarios` if intended): ",
                paste(msgs, collapse = " | ")))
  } else {
    succeed()
  }
})
