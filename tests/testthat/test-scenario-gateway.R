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

  ## Three things are diffed, and the third exists because the first two are
  ## coarse.
  ##
  ## `observed` alone is not enough. It tests `finite && total > 0`, which every
  ## scenario has satisfied since the crash fixes landed, so it is pinned at
  ## 8/8 "success" and cannot move. Measured: switching the density coordinate
  ## (#590) multiplies R0 by 2.4x to 47x and moves S01 across R0 = 1, and
  ## `observed` does not change on a single scenario. A guard that survives
  ## that is not guarding the thing the gateway now reports, so `persists` --
  ## the R0 >= 1 axis #572 made the headline -- is diffed as well.
  ##
  ## `persists` is coarse in the same way, one threshold down. Six of the eight
  ## scenarios sit at R0 below 1e-9, so they can move by orders of magnitude
  ## without crossing anything. Two changes did exactly that and neither was
  ## re-blessed: the dependency migration (#633) multiplied R0 by 8x to 32x on
  ## S03-S06 and cut S02 to 3% of its value, and the storage-pool change (#619)
  ## moved five scenarios again on top. The whole table they left behind was
  ## reported as the single flag `persists S01: TRUE -> FALSE` (#639). So
  ## `offspring_production` is diffed too, at a relative tolerance, and *any*
  ## failure prints the whole table -- one flipped flag does not distinguish a
  ## rounding change from a 32x one.
  ##
  ## The tolerance is loose enough for floating-point reordering and tight
  ## enough to catch the smallest move we have wanted to see: the 0.36% from the
  ## phylloptim 0.6.0 / odelia 0.3.1 migration (#633), which went unblessed.
  ## Raise it via SCENARIO_TOL if a platform proves noisier than that; the
  ## classifications are the part that is meant to be platform-stable.
  tol <- as.numeric(Sys.getenv("SCENARIO_TOL", "1e-3"))

  ## NA-safe inequality: a value appearing or disappearing is a change, two NAs
  ## are not. Plain `!=` gives NA there, which then indexes with NA and reports
  ## a phantom "NA: NA -> NA" row.
  differs <- function(x, y) {
    (is.na(x) != is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
  }

  report <- function(col, label) {
    if (!col %in% names(b) || !col %in% names(cur)) {
      return(NULL)
    }
    changed <- differs(b[[col]], cur[[col]])
    if (!any(changed)) {
      return(NULL)
    }
    sprintf("%s %s", label,
            paste(sprintf("%s: %s -> %s", b$scenario_id[changed],
                          b[[col]][changed], cur[[col]][changed]),
                  collapse = "; "))
  }

  ## Relative move in R0. A baseline of exactly zero has no relative scale, so
  ## it is a change iff the current value is non-zero.
  from <- b$offspring_production
  to <- cur$offspring_production
  rel <- ifelse(is.na(from) | is.na(to), NA_real_,
                ifelse(from == 0, ifelse(to == 0, 0, Inf), (to - from) / from))
  moved <- differs(from, to) & (is.na(rel) | abs(rel) > tol)

  numeric_report <- if (any(moved)) {
    sprintf("offspring_production moved on %d/%d scenario(s) (tol %g)",
            sum(moved), length(moved), tol)
  }

  msgs <- c(report("observed", "observed"), report("persists", "persists"),
            numeric_report)
  if (length(msgs)) {
    tbl <- paste(c(sprintf("%-4s %14s %14s %10s", "id", "baseline",
                           "current", "rel"),
                   sprintf("%-4s %14.7g %14.7g %10.3g%s", b$scenario_id,
                           from, to, rel, ifelse(moved, " *", ""))),
                 collapse = "\n")
    fail(paste0("Scenario outcomes changed vs baseline ",
                "(re-bless via `make bless-scenarios` if intended): ",
                paste(msgs, collapse = " | "), "\n", tbl))
  } else {
    succeed()
  }
})
