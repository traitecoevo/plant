#!/usr/bin/env Rscript
##
## Run the TF24 hydraulic scenario gateway: evaluate every scenario in
## inst/scenarios/model_scenarios_hydraulic.csv against its expected outcome and
## write a provenance-stamped scorecard.
##
## Usage:
##   Rscript scripts/run_scenario_gateway.R [output.rds] [max_patch_lifetime]
##
## Environment:
##   SCENARIO_OUT      default output path (overridden by the first CLI arg)
##   SCENARIO_MPL      default max_patch_lifetime (overridden by the second arg)
##   SCENARIO_WORKERS  parallel workers (default: detectCores() - 1)
##
## Loads the already-built shared library (no recompile); run `make compile`
## first if the C++ has changed. Parallelism is fork-based (see run_scenarios),
## which is exactly why loading via load_all() here is safe: forked workers
## inherit this session's dev namespace and compiled library. A PSOCK cluster
## would not, and would run the installed (stale) package instead.

suppressMessages(pkgload::load_all(".", recompile = FALSE, quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[[1]] else
  Sys.getenv("SCENARIO_OUT", "scenario_scorecard.rds")
mpl <- if (length(args) >= 2) as.numeric(args[[2]]) else
  as.numeric(Sys.getenv("SCENARIO_MPL", "100"))
workers <- as.integer(Sys.getenv("SCENARIO_WORKERS",
                                 as.character(max(1L, parallel::detectCores() - 1L))))

message(sprintf("Running scenarios (max_patch_lifetime = %g, workers = %d) ...",
                mpl, workers))
scorecard <- run_scenarios(max_patch_lifetime = mpl, workers = workers)

meta <- attr(scorecard, "metadata")
message(sprintf("Build: %s @ %s%s",
                substr(meta$git_commit, 1, 10), meta$git_branch,
                if (isTRUE(meta$git_dirty)) " (dirty)" else ""))

print(as.data.frame(scorecard)[, c("scenario_id", "expected", "observed",
                                   "match", "outcome", "offspring_production",
                                   "persists")],
      row.names = FALSE)

cat("\nSummary:\n")
print(as.data.frame(scenario_summary(scorecard)), row.names = FALSE)

saveRDS(scorecard, out_path)
message(sprintf("Wrote scorecard to %s", normalizePath(out_path, mustWork = FALSE)))
