# Extracted from test-scenario-gateway.R:101

# test -------------------------------------------------------------------------
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
tol <- as.numeric(Sys.getenv("SCENARIO_TOL", "1e-3"))
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
