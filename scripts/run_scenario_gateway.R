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

## Report the density coordinate: the blessed baseline depends on it (switching
## it moves R0 by up to 47x and flips S01's persistence), so a scorecard that
## does not say which one produced it cannot be compared to another.
ctrl <- scenario_control()
message(sprintf("Running scenarios (max_patch_lifetime = %g, workers = %d, %s) ...",
                mpl, workers,
                if (ctrl$node_density_in_birth_date) "density in birth date"
                else "density in height"))
scorecard <- run_scenarios(ctrl = ctrl, max_patch_lifetime = mpl,
                           workers = workers)

meta <- attr(scorecard, "metadata")
message(sprintf("Build: %s @ %s%s",
                substr(meta$git_commit, 1, 10), meta$git_branch,
                if (isTRUE(meta$git_dirty)) " (dirty)" else ""))

print(as.data.frame(scorecard)[, c("scenario_id", "expected", "observed",
                                   "match", "outcome", "offspring_production",
                                   "persists")],
      row.names = FALSE)

cat("\nSummary:\n")
s <- scenario_summary(scorecard)
print(as.data.frame(s), row.names = FALSE)

## Two axes, reported separately and labelled (#572). Deliberately NOT a single
## headline "match rate": that conflates "the model ran" with "the strategy
## lives", and with the crashes fixed (#546/#552/#554) it mostly measures how
## well the CSV's crash predictions have aged. Numerical viability is the gate the
## hydraulic work targets; persistence (R0 >= 1) is the ecological result.
cat("\n")
cat(sprintf("  Numerical viability : %d / %d ran (%.0f%%), %d crashed\n",
            s$n_ran, s$n, 100 * s$viability_rate, s$n_crashed))
cat(sprintf("  Persistence (R0>=1) : %d / %d persist (%.0f%%)\n",
            s$n_persists, s$n, 100 * s$persistence_rate))
cat(sprintf(paste("  CSV agreement       : %d / %d (%.0f%%) — agreement with the",
                  "crash-era\n                        expectations,",
                  "not a quality score\n"),
            s$n_match, s$n, 100 * s$match_rate))

saveRDS(scorecard, out_path)
message(sprintf("Wrote scorecard to %s", normalizePath(out_path, mustWork = FALSE)))
