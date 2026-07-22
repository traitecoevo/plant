#!/usr/bin/env Rscript
##
## Regenerable high-accuracy TF24 reference run.
##
## Replaces the former hand-saved `results_high.RDS`, which predated the `pars`
## refactor (#410) and several TF24 interface changes and so could no longer be
## replayed or migrated. This script produces a canonical, REPRODUCIBLE TF24
## single-species SCM run (fixed traits, fixed environment, high-accuracy
## control) for manual comparison / a coarse full-run regression anchor.
##
## Usage:
##   Rscript scripts/generate_reference_run.R [output.rds]   # default results_high.RDS
##
## Loads the already-built shared library (no recompile); run `make compile`
## first if the C++ has changed. Single SCM pass on the default node schedule
## (no adaptive refinement) with the fast default Control(), so it stays quick to
## regenerate.

suppressMessages(pkgload::load_all(".", recompile = FALSE, quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[[1]] else "results_high.RDS"

## Canonical single-species TF24 resident.
p0 <- scm_base_parameters("TF24", "TF24_Env")
p0$max_patch_lifetime <- 50
p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"))

env  <- Environment("TF24")
ctrl <- control()   # fast default (single source of truth, see NEWS §3.2)

res <- run_scm(p1, env = env, ctrl = ctrl, collect = TRUE, refine_schedule = FALSE)

saveRDS(res, out_path)
message(sprintf("model=%s  steps=%d  offspring_production=%.10g",
                model_id("TF24"), nrow(res$steps), res$offspring_production))
message(sprintf("Wrote reference run to %s",
                normalizePath(out_path, mustWork = FALSE)))
